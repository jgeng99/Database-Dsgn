/**
 * @author See Contributors.txt for code contributors and overview of BadgerDB.
 *
 * @section LICENSE
 * Copyright (c) 2012 Database Group, Computer Sciences Department, University of Wisconsin-Madison.
 */

#include "btree.h"
#include "filescan.h"
#include "exceptions/bad_index_info_exception.h"
#include "exceptions/bad_opcodes_exception.h"
#include "exceptions/bad_scanrange_exception.h"
#include "exceptions/no_such_key_found_exception.h"
#include "exceptions/scan_not_initialized_exception.h"
#include "exceptions/index_scan_completed_exception.h"
#include "exceptions/file_not_found_exception.h"
#include "exceptions/file_exists_exception.h"
#include "exceptions/end_of_file_exception.h"


//#define DEBUG

namespace badgerdb
{

// -----------------------------------------------------------------------------
// BTreeIndex::BTreeIndex -- Constructor
// -----------------------------------------------------------------------------

BTreeIndex::BTreeIndex(const std::string & relationName,
		std::string & outIndexName,
		BufMgr *bufMgrIn,
		const int attrByteOffset,
		const Datatype attrType)
{
    std::ostringstream idxStr;
    idxStr << relationName << '.' << attrByteOffset;
    outIndexName = idxStr.str(); 
    const char* meta_name = relationName.c_str();

    // initialize variables
    if (attrType==INTEGER) {
      this->leafOccupancy = INTARRAYLEAFSIZE;
      this->nodeOccupancy = INTARRAYNONLEAFSIZE;
    }
    else {
      this->leafOccupancy = (attrType==DOUBLE) ? DOUBLEARRAYLEAFSIZE:STRINGARRAYLEAFSIZE;
      this->nodeOccupancy = (attrType==DOUBLE) ? DOUBLEARRAYNONLEAFSIZE:STRINGARRAYNONLEAFSIZE;
    }

    this->attributeType = attrType;
    this->bufMgr = bufMgrIn;
    this->attrByteOffset = attrByteOffset;
    this->scanExecuting = false;

    // initialize pointers for first page of blobfile
    BlobFile* check_index;
    IndexMetaInfo* blob_meta;
    try {
        // file does not exist so create a new blob file
        check_index = new BlobFile(outIndexName, true);
        this->file = (File*) check_index;

        // allocate root and header page
        Page *blob_metaPage;
        this->bufMgr->allocPage(this->file, this->headerPageNum, blob_metaPage);

        // fill meta infor
        blob_meta = (IndexMetaInfo*) blob_metaPage;
        blob_meta->attrByteOffset = attrByteOffset;
        blob_meta->attrType = attrType;
        blob_meta->rootPageNo = rootPageNum;
        strncpy(blob_meta->relationName, meta_name, sizeof(blob_meta->relationName));

        // finish with the meta data, unpin the allocated page
        this->bufMgr->unPinPage(this->file, this->headerPageNum, true);

        // assign rootpage id through page allocation
        this->bufMgr->allocPage(this->file, this->rootPageNum, this->rootP);
        blob_meta->rootPageNo = this->rootPageNum;
        this->OrigRootId = this->rootPageNum;

        if (this->attributeType==INTEGER) {
            LeafNodeInt* rootN = (LeafNodeInt*) this->rootP;
            rootN->rightSibPageNo = 0;
            for(int i = 0; i < (this->leafOccupancy); i++) {
                rootN->keyArray[i] = 0;
            }
        }
        else if (this->attributeType==DOUBLE) {
            LeafNodeDouble* rootN = (LeafNodeDouble*) this->rootP;
            rootN->rightSibPageNo = 0;
            for(int i = 0; i < (this->leafOccupancy); i++) {
                rootN->keyArray[i] = 0;
            }
        }
        else {
            LeafNodeString* rootN = (LeafNodeString*) this->rootP;
            rootN->rightSibPageNo = 0;
            // for(int i = 0; i < (this->leafOccupancy); i++) {
            //     strncpy(rootN->keyArray[i], "", sizeof(rootN->keyArray[i]));
            // }
        }

        // unpin root page
        this->bufMgr->unPinPage(this->file, this->rootPageNum, true);

        //fill the newly created Blob File using filescan
        FileScan fileScan(relationName, bufMgr);
        RecordId scanRid;
        while(1) {
            try {
                fileScan.scanNext(scanRid);
                std::string curRecord = fileScan.getRecord();
                if ((this->attributeType==INTEGER) || 
                        (this->attributeType==DOUBLE)) {
                    insertEntry((void*) (curRecord.c_str() + attrByteOffset), scanRid);
                } 
                else {
                    char key[STRINGSIZE + 1];
                    char* src = (char*)(curRecord.c_str() + attrByteOffset);
                    strncpy(key, src, sizeof(key));
                    std::string strKey = std::string(key);
                    // std::cout<< "inserted key: "<<strKey<<std::endl;
                    this->insertEntry((void*) &strKey, scanRid);
                }
            }
            catch(EndOfFileException e) {
                break;
            }
        }
    }
    //File was not found thus we create a new one
    catch(FileExistsException e) {
        // file exists
		check_index = new BlobFile(outIndexName, false);
        this->headerPageNum = this->file->getFirstPageNo();
        this->scanExecuting = false; // no need to scan
        bool flush = false;

        this->file = (File*) check_index;

        // check metapage
        Page* blob_metaPage;
        this->bufMgr->readPage(this->file, this->headerPageNum, blob_metaPage);
        IndexMetaInfo* metaPage = (IndexMetaInfo*) blob_metaPage;
        if ((relationName.compare(metaPage->relationName)!=0) ||
            (attrByteOffset!=metaPage->attrByteOffset) ||
            (attrType!=metaPage->attrType)) {
            throw BadIndexInfoException("Meta info does not match");
        }
    
        // assign tree attributes and unpin read pages
        this->rootPageNum = metaPage->rootPageNo;
        this->bufMgr->unPinPage(this->file, this->headerPageNum, flush);
    }
}

// -----------------------------------------------------------------------------
// BTreeIndex::~BTreeIndex -- destructor
// -----------------------------------------------------------------------------

BTreeIndex::~BTreeIndex()
{
    this->scanExecuting = false;
    this->bufMgr->flushFile(this->file);
    delete this->file;
}

// -----------------------------------------------------------------------------
// BTreeIndex::insertEntry
// -----------------------------------------------------------------------------

const void BTreeIndex::insertEntry(const void *key, const RecordId rid) 
{
    this->bufMgr->readPage(this->file, this->rootPageNum, this->rootP);
    if (this->attributeType==INTEGER) {
        PageKeyPair<int>* recurPair = nullptr;
        if (!(this->rootPageNum-this->OrigRootId)) {
            this->insertRecursive(this->rootP, this->rootPageNum, true, key, rid, recurPair);
        } else this->insertRecursive(this->rootP, this->rootPageNum, false, key, rid, recurPair);
    }
    else if (this->attributeType==DOUBLE) {
       PageKeyPair<double>* recurPair = nullptr;
        if (!(this->rootPageNum-this->OrigRootId)) {
            this->insertRecursive(this->rootP, this->rootPageNum, true, key, rid, recurPair);
        } else this->insertRecursive(this->rootP, this->rootPageNum, false, key, rid, recurPair);
    }
    else {
        PageKeyPair<std::string>* recurPair = nullptr;
        if (!(this->rootPageNum-this->OrigRootId)) {
            this->insertRecursiveString(this->rootP, this->rootPageNum, true, key, rid, recurPair);
        } else {
            // std::cout<<"root change"<<std::endl;
            this->insertRecursiveString(this->rootP, this->rootPageNum, false, key, rid, recurPair);
        } 
    }
}

// -----------------------------------------------------------------------------
// BTreeIndex::insertRecursive
// -----------------------------------------------------------------------------
template <typename T>
const void BTreeIndex::insertRecursive(Page *curPage, PageId curPageNum, bool leafNode, 
                const void *key, const RecordId rid, PageKeyPair<T> *&recurPair)
{
    if (this->attributeType==INTEGER) {
        if (leafNode) {
            // leaf page helper function
            this->insertIntoLeaf(curPage, curPageNum, key, rid, recurPair);
        }
        else {
            NonLeafNodeInt *curNode = (NonLeafNodeInt*) curPage;
            // find the right key to traverse
            Page *nextPage;
            PageId nextNodeNum;
            this->findPageNoInNonLeaf(curPage, nextNodeNum, key);
            this->bufMgr->readPage(this->file, nextNodeNum, nextPage);
            leafNode = (curNode->level==1);
            this->insertRecursive(nextPage, nextNodeNum, leafNode, key, rid, recurPair);
            
            if (recurPair) {
                this->insertIntoNode(curPage, curPageNum, recurPair);
            }
            else this->bufMgr->unPinPage(this->file, curPageNum, false);
        }
    }
    else if (this->attributeType==DOUBLE) {
       if (leafNode) {
            // leaf page helper function
            this->insertIntoLeaf(curPage, curPageNum, key, rid, recurPair);
        }
        else {
            NonLeafNodeDouble *curNode = (NonLeafNodeDouble*) curPage;
            // find the right key to traverse
            Page *nextPage;
            PageId nextNodeNum;
            this->findPageNoInNonLeaf(curPage, nextNodeNum, key);
            this->bufMgr->readPage(this->file, nextNodeNum, nextPage);
            leafNode = (curNode->level==1);
            this->insertRecursive(nextPage, nextNodeNum, leafNode, key, rid, recurPair);
            
            if (recurPair) {
                this->insertIntoNode(curPage, curPageNum, recurPair);
            }
            else this->bufMgr->unPinPage(this->file, curPageNum, false);
        }
    }
}

// -----------------------------------------------------------------------------
// BTreeIndex::insertRecursiveString
// -----------------------------------------------------------------------------

const void BTreeIndex::insertRecursiveString(Page *curPage, PageId curPageNum, bool leafNode, 
                const void *key, const RecordId rid, PageKeyPair<std::string> *&recurPair)
{   
    if (leafNode) {
        // leaf page helper function
        this->insertIntoLeafString(curPage, curPageNum, key, rid, recurPair);
    }
    else {
        NonLeafNodeString *curNode = (NonLeafNodeString*) curPage;
        // find the right key to traverse
        Page *nextPage;
        PageId nextNodeNum;
        this->findPageNoInNonLeaf(curPage, nextNodeNum, key);
        this->bufMgr->readPage(this->file, nextNodeNum, nextPage);
        leafNode = (curNode->level==1);
        this->insertRecursiveString(nextPage, nextNodeNum, leafNode, key, rid, recurPair);
        
        if (recurPair) {
            this->insertIntoNodeString(curPage, curPageNum, recurPair);
        }
        else this->bufMgr->unPinPage(this->file, curPageNum, false);
    }
}

// -----------------------------------------------------------------------------
// BTreeIndex::leafFull
// -----------------------------------------------------------------------------

bool BTreeIndex::leafFull(Page* curPage) {
    if (this->attributeType==INTEGER) {
        LeafNodeInt* leafInt = (LeafNodeInt*) curPage;
        for(int i=leafOccupancy-1; i>=0; i--) {
            if(!leafInt->ridArray[i].page_number) {
                return false;
            }
        }
    }
    else if (this->attributeType==DOUBLE) {
        LeafNodeDouble* leafDouble = (LeafNodeDouble*) curPage;
        for(int i=leafOccupancy-1; i>=0; i--) {
            if(!leafDouble->ridArray[i].page_number) {
                return false;
            }
        }
    }
    else {
        LeafNodeString* leafString = (LeafNodeString*) curPage;
        for(int i=leafOccupancy-1; i>=0; i--) {
            if(!leafString->ridArray[i].page_number) {
                return false;
            }
        }
    }
    return true;
}

// -----------------------------------------------------------------------------
// BTreeIndex::nodeFull
// -----------------------------------------------------------------------------

bool BTreeIndex::nodeFull(Page* curPage) {
    if (this->attributeType==INTEGER) {
        NonLeafNodeInt* nonLeafNode = (NonLeafNodeInt*) curPage;
        for(int i=nodeOccupancy; i>=0; i--) {
            if(!nonLeafNode->pageNoArray[i]) {
                return false;
            }
        }
    }
    else if (this->attributeType==DOUBLE) {
        NonLeafNodeDouble* nonLeafNode = (NonLeafNodeDouble*) curPage;
        for(int i=nodeOccupancy; i>=0; i--) {
            if(!nonLeafNode->pageNoArray[i]) {
                return false;
            }
        }
    }
    else {
        NonLeafNodeString* nonLeafNode = (NonLeafNodeString*) curPage;
        for(int i=nodeOccupancy; i>=0; i--) {
            if(!nonLeafNode->pageNoArray[i]) {
                return false;
            }
        }
    }
    return true;
}

// -----------------------------------------------------------------------------
// BTreeIndex::insertIntoLeaf
// -----------------------------------------------------------------------------
template <typename T>
const void BTreeIndex::insertIntoLeaf(Page *curPage, PageId curPageNum,
      const void *key, const RecordId rid, PageKeyPair<T> *&recurPair) 
{
    if (!leafFull(curPage)) { // page is not full
        this->insertLeaf(curPage, key, rid);
        this->bufMgr->unPinPage(this->file, curPageNum, true);
        recurPair = nullptr;
    }
    else { // page is full, need to split
        this->splitLeaf(curPage, curPageNum, key, rid, recurPair);
    }
}

// -----------------------------------------------------------------------------
// BTreeIndex::insertIntoLeaf
// -----------------------------------------------------------------------------
const void BTreeIndex::insertIntoLeafString(Page *curPage, PageId curPageNum,
      const void *key, const RecordId rid, PageKeyPair<std::string> *&recurPair) 
{
    if (!leafFull(curPage)) { // page is not full

        this->insertLeaf(curPage, key, rid);
        this->bufMgr->unPinPage(this->file, curPageNum, true);
        recurPair = nullptr;
    }
    else { // page is full, need to split
        this->splitLeafString(curPage, curPageNum, key, rid, recurPair);
    }
}

// -----------------------------------------------------------------------------
// BTreeIndex::insertIntoNode
// -----------------------------------------------------------------------------
template <typename T>
const void BTreeIndex::insertIntoNode(Page *curPage, PageId curPageNum, PageKeyPair<T> *&recurPair) 
{
    // if page not full
    if (!nodeFull(curPage)) {
        // insert recurPair to it
        this->insertNonLeaf(curPage, recurPair);
        this->bufMgr->unPinPage(this->file, curPageNum, true);
        recurPair = nullptr;
    }
    else {
        this->splitNonLeaf(curPage, curPageNum, recurPair);
    }
}

// -----------------------------------------------------------------------------
// BTreeIndex::insertIntoNodeString
// -----------------------------------------------------------------------------

const void BTreeIndex::insertIntoNodeString(Page *curPage, PageId curPageNum, PageKeyPair<std::string> *&recurPair)
{
    // if page not full
    if (!nodeFull(curPage)) {
        // insert recurPair to it
        this->insertNonLeafString(curPage, recurPair);
        this->bufMgr->unPinPage(this->file, curPageNum, true);
        recurPair = nullptr;
    }
    else {
        this->splitNonLeafString(curPage, curPageNum, recurPair);
    }
}

// -----------------------------------------------------------------------------
// BTreeIndex::insertLeaf
// -----------------------------------------------------------------------------

const void BTreeIndex::insertLeaf(Page* curPage, const void *key, const RecordId rid)
{
    if (this->attributeType==INTEGER) {
        // determine the number of valid keys in the leaf
        LeafNodeInt* leafInt = (LeafNodeInt*) curPage;
        int keyInt = *(int*) key;
        int validKeys = 0;
        while((validKeys<leafOccupancy) && 
            (leafInt->ridArray[validKeys].page_number)) {
            validKeys++;
        }

        // find insertion point
        int start = 0;
        int end = validKeys-1;
        while(start<=end) {
            // use binary search to save runtime because sorted
            int mid = start+(end-start)/2;
            if(leafInt->keyArray[mid]<keyInt) {
                start = mid+1;
            } else {
                end = mid-1;
            }
        }
        int insertionPoint = start;

        // shift all larger keys one position to the right
        for(int i=validKeys; i>insertionPoint; i--) {
            leafInt->keyArray[i] = leafInt->keyArray[i-1];
            leafInt->ridArray[i] = leafInt->ridArray[i-1];
        }

        // insert the new key and record ID
        leafInt->keyArray[insertionPoint] = keyInt;
        leafInt->ridArray[insertionPoint] = rid;
    }
    else if (this->attributeType==DOUBLE) {
        // determine the number of valid keys in the leaf
        LeafNodeDouble* leafDouble = (LeafNodeDouble*) curPage;
        double keyDouble = *(double*) key;
        int validKeys = 0;
        while((validKeys<leafOccupancy) && 
            (leafDouble->ridArray[validKeys].page_number)) {
            validKeys++;
        }

        // find insertion point
        int start = 0;
        int end = validKeys-1;
        while(start<=end) {
            // use binary search to save runtime because sorted
            int mid = start+(end-start)/2;
            if(leafDouble->keyArray[mid]<keyDouble) {
                start = mid+1;
            } else {
                end = mid-1;
            }
        }
        int insertionPoint = start;

        // shift all larger keys one position to the right
        for(int i=validKeys; i>insertionPoint; i--) {
            leafDouble->keyArray[i] = leafDouble->keyArray[i-1];
            leafDouble->ridArray[i] = leafDouble->ridArray[i-1];
        }

        // insert the new key and record ID
        leafDouble->keyArray[insertionPoint] = keyDouble;
        leafDouble->ridArray[insertionPoint] = rid;
    }
    else {
        // determine the number of valid keys in the leaf
        LeafNodeString* leafString = (LeafNodeString*) curPage;
        std::string keyString = *((std::string*) key);
        // std::string(static_cast<const char*>(key));
        int validKeys = 0;
        while((validKeys<leafOccupancy) && 
            (leafString->ridArray[validKeys].page_number)) {
            validKeys++;
        }

        // find insertion point
        // std::cout<<"Before: "<<leafString->keyArray[0]<<std::endl<<std::endl;
        int start = 0;
        int end = validKeys-1;
        while(start<=end) {
            // use binary search to save runtime because sorted
            int mid = start+(end-start)/2;
            if(leafString->keyArray[mid]<keyString) {
                start = mid+1;
            } else {
                end = mid-1;
            }
        }
        int insertionPoint = start;

        // shift all larger keys one position to the right
        for(int i=validKeys; i>insertionPoint; i--) {
            strncpy(leafString->keyArray[i], leafString->keyArray[i-1], STRINGSIZE);
            leafString->ridArray[i] = leafString->ridArray[i-1];
        }

        // insert the new key and record ID
        strncpy(leafString->keyArray[insertionPoint], keyString.c_str(), STRINGSIZE);
        leafString->ridArray[insertionPoint] = rid;
        // std::cout<<"After: "<<leafString->keyArray[0]<<std::endl<<std::endl;
    }
}

// -----------------------------------------------------------------------------
// BTreeIndex::insertNonLeaf
// -----------------------------------------------------------------------------
template <typename T>
const void BTreeIndex::insertNonLeaf(Page* curPage, PageKeyPair<T> *recurPair)
{
    if (this->attributeType==INTEGER) {
        NonLeafNodeInt* nonLeaf = (NonLeafNodeInt*) curPage;

        // determine the number of valid keys in the non-leaf node
        int validKeys = 0;
        while((validKeys<nodeOccupancy) && 
            (nonLeaf->pageNoArray[validKeys+1])) {
            validKeys++;
        }

        // find insertion point
        int start = 0;
        int end = validKeys-1;
        while(start<=end) {
            // use binary search to save runtime because sorted
            int mid = start+(end-start)/2;
            if(nonLeaf->keyArray[mid]<recurPair->key) {
                start = mid+1;
            } else {
                end = mid-1;
            }
        }
        int insertionPoint = start;

        // shift all larger keys and page numbers one position to the right
        for(int i=validKeys; i>insertionPoint; i--) {
            nonLeaf->keyArray[i] = nonLeaf->keyArray[i-1];
            nonLeaf->pageNoArray[i+1] = nonLeaf->pageNoArray[i];
        }

        // insert the new key and corresponding page number
        nonLeaf->keyArray[insertionPoint] = recurPair->key;
        nonLeaf->pageNoArray[insertionPoint+1] = recurPair->pageNo;
    }
    else if (this->attributeType==DOUBLE) {
        NonLeafNodeDouble* nonLeaf = (NonLeafNodeDouble*) curPage;

        // determine the number of valid keys in the non-leaf node
        int validKeys = 0;
        while((validKeys<nodeOccupancy) && 
            (nonLeaf->pageNoArray[validKeys+1])) {
            validKeys++;
        }

        // find insertion point
        int start = 0;
        int end = validKeys-1;
        while(start<=end) {
            // use binary search to save runtime because sorted
            int mid = start+(end-start)/2;
            if(nonLeaf->keyArray[mid]<recurPair->key) {
                start = mid+1;
            } else {
                end = mid-1;
            }
        }
        int insertionPoint = start;

        // shift all larger keys and page numbers one position to the right
        for(int i=validKeys; i>insertionPoint; i--) {
            nonLeaf->keyArray[i] = nonLeaf->keyArray[i-1];
            nonLeaf->pageNoArray[i+1] = nonLeaf->pageNoArray[i];
        }

        // insert the new key and corresponding page number
        nonLeaf->keyArray[insertionPoint] = recurPair->key;
        nonLeaf->pageNoArray[insertionPoint+1] = recurPair->pageNo;
    }
}

// -----------------------------------------------------------------------------
// BTreeIndex::insertNonLeafString
// -----------------------------------------------------------------------------

const void BTreeIndex::insertNonLeafString(Page* curPage, PageKeyPair<std::string> *recurPair)
{
    NonLeafNodeString* nonLeaf = (NonLeafNodeString*) curPage;

    // determine the number of valid keys in the non-leaf node
    int validKeys = 0;
    while((validKeys<nodeOccupancy) && 
        (nonLeaf->pageNoArray[validKeys+1])) {
        validKeys++;
    }

    // find insertion point
    int start = 0;
    int end = validKeys-1;
    while(start<=end) {
        // use binary search to save runtime because sorted
        int mid = start+(end-start)/2;
        if(nonLeaf->keyArray[mid]<recurPair->key) {
            start = mid+1;
        } else {
            end = mid-1;
        }
    }
    int insertionPoint = start;

    // shift all larger keys and page numbers one position to the right
    for(int i=validKeys; i>insertionPoint; i--) {
        strncpy(nonLeaf->keyArray[i], nonLeaf->keyArray[i-1], STRINGSIZE);
        // nonLeaf->keyArray[i] = nonLeaf->keyArray[i-1];
        nonLeaf->pageNoArray[i+1] = nonLeaf->pageNoArray[i];
    }

    // insert the new key and corresponding page number
    strncpy(nonLeaf->keyArray[insertionPoint], recurPair->key.c_str(), STRINGSIZE);
    // std::cout<<"nonLeaf: "<<nonLeaf->keyArray[0]<<std::endl;
    // for (int i=0;i<10;i++) std::cout<<"PageNo: ("<<i<<") "<<nonLeaf->pageNoArray[i]<<" ";
    // std::cout<<std::endl;
    // nonLeaf->keyArray[insertionPoint] = recurPair->key;
    nonLeaf->pageNoArray[insertionPoint+1] = recurPair->pageNo;
}

// -----------------------------------------------------------------------------
// BTreeIndex::splitLeaf
// -----------------------------------------------------------------------------
template <typename T>
const void BTreeIndex::splitLeaf(Page* curPage, PageId leafPageId, 
          const void *key, const RecordId rid, PageKeyPair<T>*& recurPair)
{
    PageId rightId;
    Page* rightPage;
    this->bufMgr->allocPage(this->file, rightId, rightPage);

    if (this->attributeType==INTEGER) {
        // implement leaf split
        this->updateMidLeaf(curPage, rightPage, rightId, key, rid, recurPair);
        LeafNodeInt* rightNode = (LeafNodeInt*) rightPage;

        // create a new root node if the current node is a leaf
        if (leafPageId==rootPageNum) {
            // create a new root 
            PageId newParentId;
            Page* newParent;
            this->bufMgr->allocPage(this->file, newParentId, newParent);
            NonLeafNodeInt* newParentPage = (NonLeafNodeInt*) newParent;

            // level is guaranteed to be one above leaf as given
            newParentPage->level = 1;

            // set its page points to the newly split leaves
            newParentPage->pageNoArray[0] = leafPageId;
            newParentPage->pageNoArray[1] = rightId;

            // set the key of parent to be the right first key by default
            newParentPage->keyArray[0] = rightNode->keyArray[0];

            // update header page
            Page *tempMeta;
            this->bufMgr->readPage(this->file, this->headerPageNum, tempMeta);
            IndexMetaInfo* metaIndex = (IndexMetaInfo*)tempMeta;
            metaIndex->rootPageNo = newParentId;
            this->rootPageNum = newParentId;
            
            // unpin unused pages
            this->bufMgr->unPinPage(this->file, this->headerPageNum, true);
            this->bufMgr->unPinPage(this->file, newParentId, true);
        }
    }
    else if (this->attributeType==DOUBLE) {
        // implement leaf split
        this->updateMidLeaf(curPage, rightPage, rightId, key, rid, recurPair);
        LeafNodeDouble* rightNode = (LeafNodeDouble*) rightPage;

        // create a new root node if the current node is a leaf
        if (leafPageId==rootPageNum) {
            // create a new root 
            PageId newParentId;
            Page* newParent;
            this->bufMgr->allocPage(this->file, newParentId, newParent);
            NonLeafNodeDouble* newParentPage = (NonLeafNodeDouble*) newParent;

            // level is guaranteed to be one above leaf as given
            newParentPage->level = 1;

            // set its page points to the newly split leaves
            newParentPage->pageNoArray[0] = leafPageId;
            newParentPage->pageNoArray[1] = rightId;

            // set the key of parent to be the right first key by default
            newParentPage->keyArray[0] = rightNode->keyArray[0];

            // update header page
            Page *tempMeta;
            this->bufMgr->readPage(this->file, this->headerPageNum, tempMeta);
            IndexMetaInfo* metaIndex = (IndexMetaInfo*)tempMeta;
            metaIndex->rootPageNo = newParentId;
            this->rootPageNum = newParentId;
            
            // unpin unused pages
            this->bufMgr->unPinPage(this->file, this->headerPageNum, true);
            this->bufMgr->unPinPage(this->file, newParentId, true);
        }
    }

    // unpin unused pages
    this->bufMgr->unPinPage(this->file, leafPageId, true);
    this->bufMgr->unPinPage(this->file, rightId, true);
}

// -----------------------------------------------------------------------------
// BTreeIndex::splitLeafString
// -----------------------------------------------------------------------------

const void BTreeIndex::splitLeafString(Page* curPage, PageId leafPageId, 
          const void *key, const RecordId rid, PageKeyPair<std::string> *&recurPair)
{
    PageId rightId;
    Page* rightPage;
    this->bufMgr->allocPage(this->file, rightId, rightPage);

    // implement leaf split
    this->updateMidLeafString(curPage, rightPage, rightId, key, rid, recurPair);
    LeafNodeString* rightNode = (LeafNodeString*) rightPage;

    // create a new root node if the current node is a leaf

    // std::string keyOps = *((std::string*) key);
    // std::cout<<"key: "<<keyOps.substr(0,10)<<" leafPageId: "<<leafPageId<<" rootPageNum: "<<rootPageNum<<std::endl;
    if (leafPageId==rootPageNum) {
        // create a new root 
        PageId newParentId;
        Page* newParent;
        this->bufMgr->allocPage(this->file, newParentId, newParent);
        NonLeafNodeString* newParentPage = (NonLeafNodeString*) newParent;

        // level is guaranteed to be one above leaf as given
        newParentPage->level = 1;

        // set its page points to the newly split leaves
        newParentPage->pageNoArray[0] = leafPageId;
        newParentPage->pageNoArray[1] = rightId;

        // std::cout<<"parent array: "<<newParentPage->keyArray[0]<<std::endl<<std::endl;

        // set the key of parent to be the right first key by default
        strncpy(newParentPage->keyArray[0], rightNode->keyArray[0], STRINGSIZE);
        // newParentPage->keyArray[0] = rightNode->keyArray[0];

        // update header page
        Page *tempMeta;
        this->bufMgr->readPage(this->file, this->headerPageNum, tempMeta);
        IndexMetaInfo* metaIndex = (IndexMetaInfo*)tempMeta;
        metaIndex->rootPageNo = newParentId;
        this->rootPageNum = newParentId;
        
        // unpin unused pages
        this->bufMgr->unPinPage(this->file, this->headerPageNum, true);
        this->bufMgr->unPinPage(this->file, newParentId, true);
    }

    // unpin unused pages
    this->bufMgr->unPinPage(this->file, leafPageId, true);
    this->bufMgr->unPinPage(this->file, rightId, true);
}

// -----------------------------------------------------------------------------
// BTreeIndex::updateMidLeaf
// -----------------------------------------------------------------------------
template <typename T>
const void BTreeIndex::updateMidLeaf(Page* leftPage, Page* rightPage, PageId rightId,
      const void *key, const RecordId rid, PageKeyPair<T> *&recurPair)
{
    if (this->attributeType==INTEGER) {
        LeafNodeInt* leftNode = (LeafNodeInt*) leftPage;
        LeafNodeInt* rightNode = (LeafNodeInt*) rightPage;
        int keyOps = *(int*) key;

        // determine the index of the middle key
        int mid = leafOccupancy/2;

        // move the latter half of the keys and rids to the new leaf
        for (int i=0; i<leafOccupancy-mid; i++) {
            rightNode->keyArray[i] = leftNode->keyArray[i+mid];
            rightNode->ridArray[i] = leftNode->ridArray[i+mid];

            // erase the moved entries from left node
            leftNode->keyArray[i+mid] = 0;
            leftNode->ridArray[i+mid].page_number = 0;
        }

        // insert and copy up
        if (keyOps<rightNode->keyArray[0]) {
            this->insertLeaf(leftPage, key, rid);
        }
        else this->insertLeaf(rightPage, key, rid);

        // update the sibling pointers
        rightNode->rightSibPageNo = leftNode->rightSibPageNo;
        leftNode->rightSibPageNo = rightId;

        // update recurPair for the parent level
        recurPair = new PageKeyPair<T>();
        recurPair->set(rightId, rightNode->keyArray[0]);
    }
    else if (this->attributeType==DOUBLE) {
        LeafNodeDouble* leftNode = (LeafNodeDouble*) leftPage;
        LeafNodeDouble* rightNode = (LeafNodeDouble*) rightPage;
        double keyOps = *(double*) key;

        // determine the index of the middle key
        int mid = leafOccupancy/2;

        // move the latter half of the keys and rids to the new leaf
        for (int i=0; i<leafOccupancy-mid; i++) {
            rightNode->keyArray[i] = leftNode->keyArray[i+mid];
            rightNode->ridArray[i] = leftNode->ridArray[i+mid];

            // erase the moved entries from left node
            leftNode->keyArray[i+mid] = 0;
            leftNode->ridArray[i+mid].page_number = 0;
        }

        // insert and copy up
        if (keyOps<rightNode->keyArray[0]) {
            this->insertLeaf(leftPage, key, rid);
        }
        else this->insertLeaf(rightPage, key, rid);

        // update the sibling pointers
        rightNode->rightSibPageNo = leftNode->rightSibPageNo;
        leftNode->rightSibPageNo = rightId;

        // update recurPair for the parent level
        recurPair = new PageKeyPair<T>();
        recurPair->set(rightId, rightNode->keyArray[0]);
    }
}

// -----------------------------------------------------------------------------
// BTreeIndex::updateMidLeafString
// -----------------------------------------------------------------------------

const void BTreeIndex::updateMidLeafString(Page* leftPage, Page* rightPage, PageId rightId,
      const void *key, const RecordId rid, PageKeyPair<std::string> *&recurPair)
{
    LeafNodeString* leftNode = (LeafNodeString*) leftPage;
    LeafNodeString* rightNode = (LeafNodeString*) rightPage;
    std::string keyOps = *((std::string*) key);
    // std::string(static_cast<const char*>(key));

    // determine the index of the middle key
    int mid = leafOccupancy/2;


    // move the latter half of the keys and rids to the new leaf
    // std::cout<<"size of array " << sizeof(leftNode->keyArray) <<std::endl<<std::endl;
    // std::cout<<"leftArray "<<leftNode->keyArray[0]<<std::endl<<std::endl;
    for (int i=0; i<leafOccupancy-mid; i++) {
        // std::cout<<"leftArray "<<leftNode->keyArray[i]<<std::endl<<std::endl;
        strncpy(rightNode->keyArray[i], leftNode->keyArray[i+mid], STRINGSIZE);
        // rightNode->keyArray[i] = leftNode->keyArray[i+mid];
        rightNode->ridArray[i] = leftNode->ridArray[i+mid];

        // erase the moved entries from left node
        strncpy(leftNode->keyArray[i+mid], "", STRINGSIZE);
        // leftNode->keyArray[i+mid] = 0;
        leftNode->ridArray[i+mid].page_number = 0;
    }

    // insert and copy up
    // std::cout<<"leftArray "<<leftNode->keyArray[0]<<std::endl<<std::endl;
    // std::cout<<"keyarray "<<rightNode->keyArray[0]<<std::endl<<std::endl;
    // throw;
    if (keyOps<rightNode->keyArray[0]) {
        this->insertLeaf(leftPage, key, rid);
    }
    else this->insertLeaf(rightPage, key, rid);

    // update the sibling pointers
    rightNode->rightSibPageNo = leftNode->rightSibPageNo;
    leftNode->rightSibPageNo = rightId;

    // update recurPair for the parent level
    recurPair = new PageKeyPair<std::string>();
    recurPair->set(rightId, rightNode->keyArray[0]);
}

// -----------------------------------------------------------------------------
// BTreeIndex::splitNonLeaf
// -----------------------------------------------------------------------------
template <typename T>
const void BTreeIndex::splitNonLeaf(Page* curPage, PageId leafPageId, PageKeyPair<T> *&recurPair)
{
    PageId rightId;
    Page* rightPage;
    this->bufMgr->allocPage(this->file, rightId, rightPage);

    if (this->attributeType==INTEGER) {
        this->updateMidNode(curPage, rightPage, rightId, recurPair);

        // if the curNode is the root
        if (leafPageId==this->rootPageNum) {
            // create a new root 
            PageId newParentId;
            Page* newParent;
            this->bufMgr->allocPage(this->file, newParentId, newParent);
            NonLeafNodeInt* newParentPage = (NonLeafNodeInt*) newParent;

            // level can set arbitrarily as long as it != 1
            newParentPage->level = 999;

            // set its page points to the newly split leaves
            newParentPage->pageNoArray[0] = leafPageId;
            newParentPage->pageNoArray[1] = recurPair->pageNo;

            // set the key of parent to be the right first key by default
            newParentPage->keyArray[0] = recurPair->key;

            // update header page
            Page *tempMeta;
            this->bufMgr->readPage(this->file, this->headerPageNum, tempMeta);
            IndexMetaInfo* metaIndex = (IndexMetaInfo*) tempMeta;
            metaIndex->rootPageNo = newParentId;
            this->rootPageNum = newParentId;

            // unpin unused pages
            this->bufMgr->unPinPage(this->file, this->headerPageNum, true);
            this->bufMgr->unPinPage(this->file, newParentId, true);
        }
    }
    else if (this->attributeType==DOUBLE) {
        this->updateMidNode(curPage, rightPage, rightId, recurPair);

        // if the curNode is the root
        if (leafPageId==this->rootPageNum) {
            // create a new root 
            PageId newParentId;
            Page* newParent;
            this->bufMgr->allocPage(this->file, newParentId, newParent);
            NonLeafNodeDouble* newParentPage = (NonLeafNodeDouble*) newParent;

            // level can set arbitrarily as long as it != 1
            newParentPage->level = 999;

            // set its page points to the newly split leaves
            newParentPage->pageNoArray[0] = leafPageId;
            newParentPage->pageNoArray[1] = recurPair->pageNo;

            // set the key of parent to be the right first key by default
            newParentPage->keyArray[0] = recurPair->key;

            // update header page
            Page *tempMeta;
            this->bufMgr->readPage(this->file, this->headerPageNum, tempMeta);
            IndexMetaInfo* metaIndex = (IndexMetaInfo*) tempMeta;
            metaIndex->rootPageNo = newParentId;
            this->rootPageNum = newParentId;

            // unpin unused pages
            this->bufMgr->unPinPage(this->file, this->headerPageNum, true);
            this->bufMgr->unPinPage(this->file, newParentId, true);
        }
    }

    this->bufMgr->unPinPage(this->file, leafPageId, true);
    this->bufMgr->unPinPage(this->file, rightId, true);
}

// -----------------------------------------------------------------------------
// BTreeIndex::splitNonLeafString
// -----------------------------------------------------------------------------

const void BTreeIndex::splitNonLeafString(Page* curPage, PageId leafPageId, PageKeyPair<std::string> *&recurPair)
{
    PageId rightId;
    Page* rightPage;
    this->bufMgr->allocPage(this->file, rightId, rightPage);

    this->updateMidNodeString(curPage, rightPage, rightId, recurPair);

    // if the curNode is the root
    if (leafPageId==this->rootPageNum) {
        // create a new root 
        PageId newParentId;
        Page* newParent;
        this->bufMgr->allocPage(this->file, newParentId, newParent);
        NonLeafNodeString* newParentPage = (NonLeafNodeString*) newParent;

        // level can set arbitrarily as long as it != 1
        newParentPage->level = 999;

        // set its page points to the newly split leaves
        newParentPage->pageNoArray[0] = leafPageId;
        newParentPage->pageNoArray[1] = recurPair->pageNo;

        // set the key of parent to be the right first key by default
        strncpy(newParentPage->keyArray[0], recurPair->key.c_str(), STRINGSIZE);
        // newParentPage->keyArray[0] = recurPair->key;

        // update header page
        Page *tempMeta;
        this->bufMgr->readPage(this->file, this->headerPageNum, tempMeta);
        IndexMetaInfo* metaIndex = (IndexMetaInfo*) tempMeta;
        metaIndex->rootPageNo = newParentId;
        this->rootPageNum = newParentId;

        // unpin unused pages
        this->bufMgr->unPinPage(this->file, this->headerPageNum, true);
        this->bufMgr->unPinPage(this->file, newParentId, true);
    }

    this->bufMgr->unPinPage(this->file, leafPageId, true);
    this->bufMgr->unPinPage(this->file, rightId, true);
}

// -----------------------------------------------------------------------------
// BTreeIndex::updateMidNode
// -----------------------------------------------------------------------------
template <typename T>
const void BTreeIndex::updateMidNode(Page* leftPage, Page* rightPage, PageId rightId,
                  PageKeyPair<T> *&recurPair)
{
    // std::cout << "------------updateMidLeaf used------------" << std::endl;
    if (this->attributeType==INTEGER) {
        NonLeafNodeInt* leftNode = (NonLeafNodeInt*) leftPage;
        NonLeafNodeInt* rightNode = (NonLeafNodeInt*) rightPage;

        // determine the index of the middle key, increment by one because push up
        int mid = leafOccupancy/2 + 1;

        // move the latter half of the keys and rids to the new leaf
        for(int i=0; i<nodeOccupancy-mid; i++) {
            rightNode->keyArray[i] = leftNode->keyArray[i+mid];
            rightNode->pageNoArray[i] = leftNode->pageNoArray[i+mid+1];

            // erase the moved entries from left node
            leftNode->keyArray[i+mid] = 0;
            leftNode->pageNoArray[i+mid+1] = 0;
        }

        // update recurPair for the parent level
        recurPair = new PageKeyPair<T>();
        recurPair->set(rightId, leftNode->keyArray[mid-1]);

        // push up instead of copy up
        leftNode->keyArray[mid-1] = 0;
        leftNode->pageNoArray[mid] = 0;

        // insert the new child entry
        if (recurPair->key<rightNode->keyArray[0]) {
            this->insertNonLeaf(leftPage, recurPair);
        } 
        else this->insertNonLeaf(rightPage, recurPair);

        // update the level attribute
        rightNode->level = leftNode->level;
    }
    else if (this->attributeType==DOUBLE) {
        NonLeafNodeDouble* leftNode = (NonLeafNodeDouble*) leftPage;
        NonLeafNodeDouble* rightNode = (NonLeafNodeDouble*) rightPage;

        // determine the index of the middle key, increment by one because push up
        int mid = leafOccupancy/2 + 1;

        // move the latter half of the keys and rids to the new leaf
        for(int i=0; i<nodeOccupancy-mid; i++) {
            rightNode->keyArray[i] = leftNode->keyArray[i+mid];
            rightNode->pageNoArray[i] = leftNode->pageNoArray[i+mid+1];

            // erase the moved entries from left node
            leftNode->keyArray[i+mid] = 0;
            leftNode->pageNoArray[i+mid+1] = 0;
        }

        // update recurPair for the parent level
        recurPair = new PageKeyPair<T>();
        recurPair->set(rightId, leftNode->keyArray[mid-1]);

        // push up instead of copy up
        leftNode->keyArray[mid-1] = 0;
        leftNode->pageNoArray[mid] = 0;

        // insert the new child entry
        if (recurPair->key<rightNode->keyArray[0]) {
            this->insertNonLeaf(leftPage, recurPair);
        } 
        else this->insertNonLeaf(rightPage, recurPair);

        // update the level attribute
        rightNode->level = leftNode->level;
    }
}

// -----------------------------------------------------------------------------
// BTreeIndex::updateMidNodeString
// -----------------------------------------------------------------------------
const void BTreeIndex::updateMidNodeString(Page* leftPage, Page* rightPage, PageId rightId,
                  PageKeyPair<std::string> *&recurPair)
{
    NonLeafNodeString* leftNode = (NonLeafNodeString*) leftPage;
    NonLeafNodeString* rightNode = (NonLeafNodeString*) rightPage;

    // determine the index of the middle key, increment by one because push up
    int mid = leafOccupancy/2 + 1;

    // move the latter half of the keys and rids to the new leaf
    for(int i=0; i<nodeOccupancy-mid; i++) {
        strncpy(rightNode->keyArray[i], leftNode->keyArray[i+mid], STRINGSIZE);
        // rightNode->keyArray[i] = leftNode->keyArray[i+mid];
        rightNode->pageNoArray[i] = leftNode->pageNoArray[i+mid+1];

        // erase the moved entries from left node
        strncpy(leftNode->keyArray[i+mid], "", STRINGSIZE);
        // leftNode->keyArray[i+mid] = 0;
        leftNode->pageNoArray[i+mid+1] = 0;
    }

    // update recurPair for the parent level
    recurPair = new PageKeyPair<std::string>();
    recurPair->set(rightId, leftNode->keyArray[mid-1]);

    // push up instead of copy up
    strncpy(leftNode->keyArray[mid-1], "", STRINGSIZE);
    // leftNode->keyArray[mid-1] = 0;
    leftNode->pageNoArray[mid] = 0;

    // insert the new child entry
    if (recurPair->key<rightNode->keyArray[0]) {
        this->insertNonLeafString(leftPage, recurPair);
    } 
    else this->insertNonLeafString(rightPage, recurPair);

    // update the level attribute
    rightNode->level = leftNode->level;
}

// -----------------------------------------------------------------------------
// BTreeIndex::findPageNoInNonLeaf
// -----------------------------------------------------------------------------

const void BTreeIndex::findPageNoInNonLeaf(Page *curPage, PageId &nextNodeNum, const void *key)
{
    if (this->attributeType==INTEGER) {
        NonLeafNodeInt* nonLeafNode = (NonLeafNodeInt*) curPage;
        int keyOp = *(int*) key;
        for (int i=nodeOccupancy-1; i>=0; i--) {
            if (nonLeafNode->pageNoArray[i]) {
                if (!i || nonLeafNode->keyArray[i-1]<keyOp) {
                    nextNodeNum = nonLeafNode->pageNoArray[i];
                    return;
                }
            }
        }
    }
    else if (this->attributeType==DOUBLE) {
        NonLeafNodeDouble* nonLeafNode = (NonLeafNodeDouble*) curPage;
        double keyOp = *(double*) key;
        for (int i=nodeOccupancy-1; i>=0; i--) {
            if (nonLeafNode->pageNoArray[i]) {
                if (!i || nonLeafNode->keyArray[i-1]<keyOp) {
                    nextNodeNum = nonLeafNode->pageNoArray[i];
                    return;
                }
            }
        }
    }
    else {
        std::string keyOp = std::string(static_cast<const char*>(key)); 
        NonLeafNodeString* nonLeafNode = (NonLeafNodeString*) curPage;
        // std::string keyOp = *(std::string*) key;

        // char keyChar[STRINGSIZE + 1];
        // char* src = (char*)(key + this->attrByteOffset);
        // strncpy(keyChar, src, sizeof(keyChar));
        // std::string keyOp = std::string(keyChar);

        for (int i=nodeOccupancy-1; i>=0; i--) {
            if (nonLeafNode->pageNoArray[i]) {
                // std::string keyOp = *(std::string*) key;
                if (!i || nonLeafNode->keyArray[i-1]<keyOp) {
                    std::cout<<"key: "<<keyOp<<" Go to page: "<<i<<std::endl;
                    nextNodeNum = nonLeafNode->pageNoArray[i];
                    return;
                }
            }
        }
    }
}

// -----------------------------------------------------------------------------
// BTreeIndex::startScan
// -----------------------------------------------------------------------------

const void BTreeIndex::startScan(const void* lowValParm,
           const Operator lowOpParm,
           const void* highValParm,
           const Operator highOpParm)
{
    if (this->scanExecuting) endScan();
    if (lowOpParm!=GT && lowOpParm!=GTE) throw BadOpcodesException();
    if (highOpParm!=LT && highOpParm!=LTE) throw BadOpcodesException();

    if (this->attributeType==INTEGER) {
        // sanity check
        this->lowValInt = *(int*)lowValParm;
        this->highValInt = *(int*)highValParm;
        if (lowValInt>highValInt) throw BadScanrangeException();
        this->lowOp = lowOpParm;
        this->highOp = highOpParm;
        this->currentPageNum = this->rootPageNum;
        this->bufMgr->readPage(this->file, this->currentPageNum, this->currentPageData);

        // root is not a leaf
        if(this->rootPageNum-this->OrigRootId) {
            this->findStartLeaf(lowValParm);
        }
        // find the smallest one that meet the predicate on leaf
        bool keyFound = false;
        this->findKeyLeaf(keyFound);
        // std::cout <<"--------"<<found<<"--------"<<std::endl;
        if (!keyFound) {
            this->bufMgr->unPinPage(this->file, this->currentPageNum, false);
            throw NoSuchKeyFoundException();
        }
    }
    else if (this->attributeType==DOUBLE) {
        // sanity check
        this->lowValDouble = *(double*)lowValParm;
        this->highValDouble = *(double*)highValParm;
        if (lowValDouble>highValDouble) throw BadScanrangeException();
        this->lowOp = lowOpParm;
        this->highOp = highOpParm;
        this->currentPageNum = this->rootPageNum;
        this->bufMgr->readPage(this->file, this->currentPageNum, this->currentPageData);

        // root is not a leaf
        if(this->rootPageNum-this->OrigRootId) {
            this->findStartLeaf(lowValParm);
        }
        // find the smallest one that meet the predicate on leaf
        bool keyFound = false;
        this->findKeyLeaf(keyFound);
        // std::cout <<"--------"<<found<<"--------"<<std::endl;
        if (!keyFound) {
            this->bufMgr->unPinPage(this->file, this->currentPageNum, false);
            throw NoSuchKeyFoundException();
        }
    }
    else {
        // sanity check
        this->lowValString = std::string(static_cast<const char*>(lowValParm));
        this->highValString = std::string(static_cast<const char*>(highValParm));
        this->lowValString = this->lowValString.substr(0,STRINGSIZE);
        this->highValString = this->highValString.substr(0,STRINGSIZE);

        if (lowValString>highValString) throw BadScanrangeException();
        this->lowOp = lowOpParm;
        this->highOp = highOpParm;
        this->currentPageNum = this->rootPageNum;
        this->bufMgr->readPage(this->file, this->currentPageNum, this->currentPageData);

        // root is not a leaf
        // std::cout<<this->rootPageNum<<std::endl<<std::endl;
        // std::cout<<this->OrigRootId<<std::endl<<std::endl;
        if(this->rootPageNum-this->OrigRootId) {
            // std::cout<<"12345"<<std::endl<<std::endl;
            this->findStartLeaf(lowValParm);
        }
        // find the smallest one that meet the predicate on leaf
        bool keyFound = false;
        this->findKeyLeaf(keyFound);
        // std::cout <<"--------"<<found<<"--------"<<std::endl;
        if (!keyFound) {
            this->bufMgr->unPinPage(this->file, this->currentPageNum, false);
            throw NoSuchKeyFoundException();
        }
    }
}

// -----------------------------------------------------------------------------
// BTreeIndex::findStartLeaf
// -----------------------------------------------------------------------------

const void BTreeIndex::findStartLeaf(const void* lowValParm) {
    if (this->attributeType==INTEGER) {
        NonLeafNodeInt* currentNode = (NonLeafNodeInt*) this->currentPageData;
        PageId nextPageNum;

        // reach toward the parent of a leaf
        while(currentNode->level!=1) {
            // find the child page given some predicate
            this->findPageNoInNonLeaf(this->currentPageData, nextPageNum, lowValParm);
            this->bufMgr->unPinPage(this->file, this->currentPageNum, false);

            // proceed to the child page
            this->currentPageNum = nextPageNum;
            this->bufMgr->readPage(this->file, this->currentPageNum, this->currentPageData);
            currentNode = (NonLeafNodeInt*) this->currentPageData;
        }

        // parent of a leaf
        this->findPageNoInNonLeaf(this->currentPageData, nextPageNum, lowValParm);
        this->bufMgr->unPinPage(this->file, this->currentPageNum, false);
        this->currentPageNum = nextPageNum;
        this->bufMgr->readPage(this->file, this->currentPageNum, this->currentPageData);
    }
    else if (this->attributeType==DOUBLE) {
        NonLeafNodeDouble* currentNode = (NonLeafNodeDouble*) this->currentPageData;
        PageId nextPageNum;

        // reach toward the parent of a leaf
        while(currentNode->level!=1) {
            // find the child page given some predicate
            this->findPageNoInNonLeaf(this->currentPageData, nextPageNum, lowValParm);
            this->bufMgr->unPinPage(this->file, this->currentPageNum, false);

            // proceed to the child page
            this->currentPageNum = nextPageNum;
            this->bufMgr->readPage(this->file, this->currentPageNum, this->currentPageData);
            currentNode = (NonLeafNodeDouble*) this->currentPageData;
        }

        // parent of a leaf
        this->findPageNoInNonLeaf(this->currentPageData, nextPageNum, lowValParm);
        this->bufMgr->unPinPage(this->file, this->currentPageNum, false);
        this->currentPageNum = nextPageNum;
        this->bufMgr->readPage(this->file, this->currentPageNum, this->currentPageData);
    }
    else {
        NonLeafNodeString* currentNode = (NonLeafNodeString*) this->currentPageData;
        PageId nextPageNum;

        // reach toward the parent of a leaf
        while(currentNode->level!=1) {
            // find the child page given some predicate
            this->findPageNoInNonLeaf(this->currentPageData, nextPageNum, lowValParm);
            this->bufMgr->unPinPage(this->file, this->currentPageNum, false);

            // proceed to the child page
            this->currentPageNum = nextPageNum;
            this->bufMgr->readPage(this->file, this->currentPageNum, this->currentPageData);
            currentNode = (NonLeafNodeString*) this->currentPageData;
        }

        // parent of a leaf
        // std::cout<<"Before: "<<this->currentPageNum<<std::endl;
        this->findPageNoInNonLeaf(this->currentPageData, nextPageNum, lowValParm);
        this->bufMgr->unPinPage(this->file, this->currentPageNum, false);
        this->currentPageNum = nextPageNum;
        // std::cout<<"After: "<<this->currentPageNum<<std::endl;
        this->bufMgr->readPage(this->file, this->currentPageNum, this->currentPageData);
    }
}

// -----------------------------------------------------------------------------
// BTreeIndex::findKeyLeaf
// -----------------------------------------------------------------------------

const void BTreeIndex::findKeyLeaf(bool& keyFound) {
    // find the smallest one that meet the predicate
    if (this->attributeType==INTEGER) {
        LeafNodeInt* currentNode = (LeafNodeInt*) this->currentPageData;

        while ((currentNode->ridArray[0].page_number)) {
            // search key on one page
            for(int i = 0; i < leafOccupancy; i++) {
                int key = currentNode->keyArray[i];

                // if key has been inserted
                if(!currentNode->ridArray[i].page_number) {
                    break;
                }

                // find key
                if (((lowOp==GT && key>lowValInt) || (lowOp==GTE && key>=lowValInt)) && 
                            ((highOp==LT && key<highValInt) || (highOp==LTE && key<=highValInt))) {
                    this->nextEntry = i;
                    this->scanExecuting = true;
                    keyFound = true;
                    return;
                }

                // check edge cases
                else if ((highOp==LT && key>=highValInt) 
                        || (highOp==LTE && key>highValInt)) {
                    return;
                }
            }
            // no next leaf
            if(!currentNode->rightSibPageNo) {
                return;
            }

            // move to the next leaf
            this->bufMgr->unPinPage(this->file, this->currentPageNum, false);
            this->currentPageNum = currentNode->rightSibPageNo;
            this->bufMgr->readPage(this->file, this->currentPageNum, this->currentPageData);
            currentNode = (LeafNodeInt*) this->currentPageData;
        }
    }
    else if (this->attributeType==DOUBLE) {
        LeafNodeDouble* currentNode = (LeafNodeDouble*) this->currentPageData;

        while ((currentNode->ridArray[0].page_number)) {
            // search key on one page
            for(int i = 0; i < leafOccupancy; i++) {
                double key = currentNode->keyArray[i];

                // if key has been inserted
                if(!currentNode->ridArray[i].page_number) {
                    break;
                }

                // find key
                if (((lowOp==GT && key>lowValDouble) || (lowOp==GTE && key>=lowValDouble)) && 
                            ((highOp==LT && key<highValDouble) || (highOp==LTE && key<=highValDouble))) {
                    this->nextEntry = i;
                    this->scanExecuting = true;
                    keyFound = true;
                    return;
                }

                // check edge cases
                else if ((highOp==LT && key>=highValDouble) 
                        || (highOp==LTE && key>highValDouble)) {
                    return;
                }
            }
            // no next leaf
            if(!currentNode->rightSibPageNo) {
                return;
            }

            // move to the next leaf
            this->bufMgr->unPinPage(this->file, this->currentPageNum, false);
            this->currentPageNum = currentNode->rightSibPageNo;
            this->bufMgr->readPage(this->file, this->currentPageNum, this->currentPageData);
            currentNode = (LeafNodeDouble*) this->currentPageData;
        }
    }
    else {
        LeafNodeString* currentNode = (LeafNodeString*) this->currentPageData;

        // std::cout<<"starting keys: "<<currentNode->keyArray[0]<<std::endl<<std::endl;
        while ((currentNode->ridArray[0].page_number)) {
            // search key on one page
            for(int i = 0; i < leafOccupancy; i++) {
                // no page
                if(!currentNode->ridArray[i].page_number) {
                    break;
                }
                // std::string(static_cast<const char*>(lowValParm))
                std::string key = currentNode->keyArray[i];
                key = key.substr(0, STRINGSIZE);
                // std::cout<<this->nextEntry<<"lowValString("<<lowValString<<") < key("<<key<<") < highstring(";
                // std::cout<<highValString<<") "<<(key<highValString&&key>lowValString)<<std::endl;

                // find key
                if (((lowOp==GT && key>lowValString) 
                    || (lowOp==GTE && key>=lowValString)) 
                    && ((highOp==LT && key<highValString) 
                    || (highOp==LTE && key<=highValString))) {
                    this->nextEntry = i;
                    // std::cout<<currentNode->keyArray[this->nextEntry]<<std::endl;
                    this->scanExecuting = true;
                    keyFound = true;
                    // std::cout<<"the starting key is: "<<key<<std::endl<<std::endl;
                    return;
                }

                // check edge cases
                else if ((highOp==LT && key>=highValString) 
                        || (highOp==LTE && key>highValString)) {
                    return;
                }
            }
            // no next leaf
            if(!currentNode->rightSibPageNo) {
                return;
            }

            // move to the next leaf
            this->bufMgr->unPinPage(this->file, this->currentPageNum, false);
            this->currentPageNum = currentNode->rightSibPageNo;
            this->bufMgr->readPage(this->file, this->currentPageNum, this->currentPageData);
            currentNode = (LeafNodeString*) this->currentPageData;
        }
    }
}

// -----------------------------------------------------------------------------
// BTreeIndex::scanNext
// -----------------------------------------------------------------------------

const void BTreeIndex::scanNext(RecordId& outRid) 
{
    if(!this->scanExecuting) throw ScanNotInitializedException();
    if (this->attributeType==INTEGER) {
        LeafNodeInt* currentNode = (LeafNodeInt*) this->currentPageData;

        if (this->nextEntry<leafOccupancy) {
            if (currentNode->ridArray[this->nextEntry].page_number) {
                int key = currentNode->keyArray[this->nextEntry];
                if (((lowOp==GT && key>lowValInt) || (lowOp==GTE && key>=lowValInt)) && 
                        ((highOp==LT && key<highValInt) || (highOp==LTE && key<=highValInt))) {
                    outRid = currentNode->ridArray[this->nextEntry];
                    this->nextEntry++;
                    return;
                }
            }
            else if (!currentNode->rightSibPageNo) throw IndexScanCompletedException();
        }

        // change to right leaf
        this->bufMgr->unPinPage(this->file, this->currentPageNum, false);
        this->currentPageNum = currentNode->rightSibPageNo;
        this->bufMgr->readPage(this->file, currentPageNum, this->currentPageData);
        currentNode = (LeafNodeInt*) this->currentPageData;
        this->nextEntry = 0;

        // check again
        int key = currentNode->keyArray[this->nextEntry];
        if (((lowOp==GT && key>lowValInt) || (lowOp==GTE && key>=lowValInt)) && 
                ((highOp==LT && key<highValInt) || (highOp==LTE && key<=highValInt))) {
            outRid = currentNode->ridArray[this->nextEntry];
            this->nextEntry++;
            return;
        }
    }
    else if (this->attributeType==DOUBLE) {
        LeafNodeDouble* currentNode = (LeafNodeDouble*) this->currentPageData;

        if (this->nextEntry<leafOccupancy) {
            if (currentNode->ridArray[this->nextEntry].page_number) {
                double key = currentNode->keyArray[this->nextEntry];
                if (((lowOp==GT && key>lowValDouble) || (lowOp==GTE && key>=lowValDouble)) && 
                        ((highOp==LT && key<highValDouble) || (highOp==LTE && key<=highValDouble))) {
                    outRid = currentNode->ridArray[this->nextEntry];
                    this->nextEntry++;
                    return;
                }
            }
            else if (!currentNode->rightSibPageNo) throw IndexScanCompletedException();
        }

        // change to right leaf
        this->bufMgr->unPinPage(this->file, this->currentPageNum, false);
        this->currentPageNum = currentNode->rightSibPageNo;
        this->bufMgr->readPage(this->file, currentPageNum, this->currentPageData);
        currentNode = (LeafNodeDouble*) this->currentPageData;
        this->nextEntry = 0;

        // check again
        double key = currentNode->keyArray[this->nextEntry];
        if (((lowOp==GT && key>lowValDouble) || (lowOp==GTE && key>=lowValDouble)) && 
                ((highOp==LT && key<highValDouble) || (highOp==LTE && key<=highValDouble))) {
            outRid = currentNode->ridArray[this->nextEntry];
            this->nextEntry++;
            return;
        }
    }
    else {
        LeafNodeString* currentNode = (LeafNodeString*) this->currentPageData;

        if (this->nextEntry<leafOccupancy) {
            // Page* temp;
            // this->bufMgr->readPage(this->file, currentNode->rightSibPageNo, temp);
            // LeafNodeString* tempN = (LeafNodeString*) temp;
            // std::cout<<tempN->keyArray[0];
            // throw;
            if (currentNode->ridArray[this->nextEntry].page_number) {
                std::string key = currentNode->keyArray[this->nextEntry];
                key = key.substr(0,STRINGSIZE);
                // std::cout<<nextEntry<<"lowValString("<<lowValString<<") < key("<<key<<") < highstring(";
                // std::cout<<highValString<<") "<<(key<highValString&&key>lowValString)<<std::endl;
                if (((lowOp==GT && key>lowValString) 
                    || (lowOp==GTE && key>=lowValString)) 
                    && ((highOp==LT && key<highValString) 
                    || (highOp==LTE && key<=highValString))) {
                    outRid = currentNode->ridArray[this->nextEntry];
                    this->nextEntry++;
                    return;
                }
            }
            else if (!currentNode->rightSibPageNo) {
                std::cout << "End of file!" << std::endl;
                throw IndexScanCompletedException();
            }
        }

        // change to right leaf
        this->bufMgr->unPinPage(this->file, this->currentPageNum, false);
        this->currentPageNum = currentNode->rightSibPageNo;
        this->bufMgr->readPage(this->file, currentPageNum, this->currentPageData);
        currentNode = (LeafNodeString*) this->currentPageData;
        this->nextEntry = 0;

        // check again
        std::string key = currentNode->keyArray[this->nextEntry];
        key = key.substr(0,STRINGSIZE);
        // std::cout<<nextEntry<<"lowValString("<<lowValString<<") < key("<<key<<") < highstring(";
        // std::cout<<highValString<<") "<<(key<highValString&&key>lowValString)<<std::endl;
        if (((lowOp==GT && key>lowValString) 
            || (lowOp==GTE && key>=lowValString)) 
            && ((highOp==LT && key<highValString) 
            || (highOp==LTE && key<=highValString))) {
            outRid = currentNode->ridArray[this->nextEntry];
            this->nextEntry++;
            return;
        }
    }

    throw IndexScanCompletedException();
}

// -----------------------------------------------------------------------------
// BTreeIndex::endScan
// -----------------------------------------------------------------------------

const void BTreeIndex::endScan() 
{
    if(!this->scanExecuting) throw ScanNotInitializedException();
    if (this->currentPageNum) this->bufMgr->unPinPage(this->file, this->currentPageNum, false);
    this->scanExecuting = false;
}

}