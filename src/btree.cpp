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
    /* VARIABLES:
      1: bufMgr, attributeType, attrByteOffset, leafOccupancy, nodeOccupancy,
        file, headerPageNum, rootPageNum
    */

    /* VARIABLES SCANNING
      1: scanExecuting, 

      2: nextEntry, currentPageNum, currentPageData,
      lowValInt (double/string), highValInt (double/string),
      lowOp, highOp
    */
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
    //File was not found thus we create a new one
    catch(FileNotFoundException e) {
        //File did not exist from upon, thus create a new blob file
        check_index = new BlobFile(outIndexName, true);
        this->file = (File*) check_index;
        this->heightTree = 0;

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
        this->initialRootPageNum = this->rootPageNum;

        if (this->attributeType==INTEGER) {
            LeafNodeInt* rootN = (LeafNodeInt*) this->rootP;
            rootN->rightSibPageNo = 0;
            for(int i = 0; i < (this->leafOccupancy); i++) {
                rootN->keyArray[i] = 0;
            }
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
                // else {
                        // 	char key[STRINGSIZE + 1];
                //             char* src = (char*)(curRecord.c_str() + attrByteOffset);
                //             strncpy(key, src, sizeof(key));
                        // 	std::string strKey = std::string(key);
                        // 	this->insertEntry((void*) &strKey, scanRid);
                        // }
            }
            catch(EndOfFileException e) {
                break;
            }
        }
    }
}

// -----------------------------------------------------------------------------
// BTreeIndex::~BTreeIndex -- destructor
// -----------------------------------------------------------------------------

BTreeIndex::~BTreeIndex()
{
    this->scanExecuting = false;
    this->bufMgr->flushFile(this->file);
    delete file;
}

// -----------------------------------------------------------------------------
// BTreeIndex::insertEntry
// -----------------------------------------------------------------------------

const void BTreeIndex::insertEntry(const void *key, const RecordId rid) 
{
  bufMgr->readPage(file, rootPageNum, rootP);
  recurPair = nullptr;

  if (initialRootPageNum == rootPageNum) {
      insertRecursive(rootP, rootPageNum, true, key, rid, recurPair);
  } else insertRecursive(rootP, rootPageNum, false, key, rid, recurPair);
}

// -----------------------------------------------------------------------------
// BTreeIndex::insertEntry
// -----------------------------------------------------------------------------

const void BTreeIndex::insertRecursive(Page *curPage, PageId curPageNum, bool nodeIsLeaf, 
                const void *key, const RecordId rid, PageKeyPair<int> *&recurPair)
{
  if (nodeIsLeaf)
  {
      // leaf page helper function
      insertIntoLeaf(curPage, curPageNum, key, rid, recurPair);
  }
  else
  {
      NonLeafNodeInt *curNode = (NonLeafNodeInt *)curPage;
      // find the right key to traverse
      Page *nextPage;
      PageId nextNodeNum;
      findPageNoInNonLeaf(curNode, nextNodeNum, key);
      bufMgr->readPage(file, nextNodeNum, nextPage);
      nodeIsLeaf = curNode->level == 1;
      insertRecursive(nextPage, nextNodeNum, nodeIsLeaf, key, rid, recurPair);
    
      if (recurPair) {
          insertIntoNode(curPage, curPageNum, recurPair);
      }
      else bufMgr->unPinPage(file, curPageNum, false);
  }
}

// -----------------------------------------------------------------------------
// BTreeIndex::leafFull
// -----------------------------------------------------------------------------

bool BTreeIndex::leafFull(Page* curPage) {
    LeafNodeInt* leafInt = (LeafNodeInt*) curPage;
    for(int i=leafOccupancy-1; i >= 0; i--) {
        if(leafInt->ridArray[i].page_number == 0) {
            return false;
        }
    }
    return true;
}

// -----------------------------------------------------------------------------
// BTreeIndex::nodeFull
// -----------------------------------------------------------------------------

bool BTreeIndex::nodeFull(Page* curPage) {
    NonLeafNodeInt* nonLeafNode = (NonLeafNodeInt*) curPage;
    for(int i=nodeOccupancy; i >= 0; i--) {
        if(nonLeafNode->pageNoArray[i] == 0) {
            return false;
        }
    }
    return true;
}

// -----------------------------------------------------------------------------
// BTreeIndex::insertIntoLeaf
// -----------------------------------------------------------------------------

const void BTreeIndex::insertIntoLeaf(Page *curPage, PageId curPageNum,
      const void *key, const RecordId rid, PageKeyPair<int> *&recurPair) 
{
    if (!leafFull(curPage))
    { // page is not full
        insertLeaf(curPage, key, rid);
        bufMgr->unPinPage(file, curPageNum, true);
        recurPair = nullptr;
    }
    else
    { // page is full, need to split
        splitLeaf(curPage, curPageNum, key, rid, recurPair);
    }
}

// -----------------------------------------------------------------------------
// BTreeIndex::insertIntoNode
// -----------------------------------------------------------------------------

const void BTreeIndex::insertIntoNode(Page *curPage, PageId curPageNum, PageKeyPair<int> *&recurPair) 
{
    // if page not full
    if (!nodeFull(curPage))
    {
        // insert recurPair to it
        insertNonLeaf(curPage, recurPair);
        recurPair = nullptr;
        bufMgr->unPinPage(file, curPageNum, true);
    }
    else
    {
        splitNonLeaf(curPage, curPageNum, recurPair);
    }
}

// -----------------------------------------------------------------------------
// BTreeIndex::insertLeaf
// -----------------------------------------------------------------------------

const void BTreeIndex::insertLeaf(Page* curPage, const void *key, const RecordId rid)
{
    // determine the number of valid keys in the leaf
    LeafNodeInt* leafInt = (LeafNodeInt*) curPage;
    int keyInt = *(int*) key;
    int validKeys = 0;
    while((validKeys<leafOccupancy) && 
        (leafInt->ridArray[validKeys].page_number!=0)) {
        validKeys++;
    }

    // find insertion point
    int start = 0;
    int end = validKeys - 1;
    while(start <= end) {
        // use binary search to save runtime because sorted
        int mid = start + (end - start) / 2;
        if(leafInt->keyArray[mid] < keyInt) {
            start = mid + 1;
        } else {
            end = mid - 1;
        }
    }
    int insertionPoint = start;

    // shift all larger keys one position to the right
    for(int i = validKeys; i > insertionPoint; i--) {
        leafInt->keyArray[i] = leafInt->keyArray[i-1];
        leafInt->ridArray[i] = leafInt->ridArray[i-1];
    }

    // insert the new key and record ID
    leafInt->keyArray[insertionPoint] = keyInt;
    leafInt->ridArray[insertionPoint] = rid;
}

// -----------------------------------------------------------------------------
// BTreeIndex::insertNonLeaf
// -----------------------------------------------------------------------------

const void BTreeIndex::insertNonLeaf(Page* curPage, PageKeyPair<int> *recurPair)
{
    NonLeafNodeInt* nonLeaf = (NonLeafNodeInt*) curPage;

    // determine the number of valid keys in the non-leaf node
    int validKeys = 0;
    while((validKeys<nodeOccupancy) && 
        (nonLeaf->pageNoArray[validKeys+1]!=0)) {
        validKeys++;
    }

    // find insertion point
    int start = 0;
    int end = validKeys - 1;
    while(start <= end) {
        // use binary search to save runtime because sorted
        int mid = start + (end - start) / 2;
        if(nonLeaf->keyArray[mid] < recurPair->key) {
            start = mid + 1;
        } else {
            end = mid - 1;
        }
    }
    int insertionPoint = start;

    // shift all larger keys and page numbers one position to the right
    for(int i = validKeys; i > insertionPoint; i--) {
        nonLeaf->keyArray[i] = nonLeaf->keyArray[i-1];
        nonLeaf->pageNoArray[i+1] = nonLeaf->pageNoArray[i];
    }

    // insert the new key and corresponding page number
    nonLeaf->keyArray[insertionPoint] = recurPair->key;
    nonLeaf->pageNoArray[insertionPoint + 1] = recurPair->pageNo;
}

// -----------------------------------------------------------------------------
// BTreeIndex::splitLeaf
// -----------------------------------------------------------------------------

const void BTreeIndex::splitLeaf(Page* curPage, PageId leafPageId, 
          const void *key, const RecordId rid, PageKeyPair<int>*& recurPair)
{
    PageId rightId;
    Page* rightPage;
    bufMgr->allocPage(file, rightId, rightPage);

    // implement leaf split
    updateMidLeaf(curPage, rightPage, rightId, key, rid, recurPair);
    LeafNodeInt* rightNode = (LeafNodeInt*) rightPage;

    // create a new root node if the current node is a leaf
    if (leafPageId == rootPageNum)
    {
        // create a new root 
        PageId newParentId;
        Page* newParent;
        bufMgr->allocPage(file, newParentId, newParent);
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

    // Unpin unused pages
    this->bufMgr->unPinPage(this->file, leafPageId, true);
    this->bufMgr->unPinPage(this->file, rightId, true);
}

// -----------------------------------------------------------------------------
// BTreeIndex::updateMidLeaf
// -----------------------------------------------------------------------------

const void BTreeIndex::updateMidLeaf(Page* leftPage, Page* rightPage, PageId rightId,
      const void *key, const RecordId rid, PageKeyPair<int> *&recurPair)
{
    LeafNodeInt* leftNode = (LeafNodeInt*) leftPage;
    LeafNodeInt* rightNode = (LeafNodeInt*) rightPage;
    int keyOps = *(int*) key;

    // determine the index of the middle key
    int mid = leafOccupancy / 2;

    // move the latter half of the keys and rids to the new leaf
    for (int i = 0; i < leafOccupancy-mid; i++) {
        rightNode->keyArray[i] = leftNode->keyArray[i+mid];
        rightNode->ridArray[i] = leftNode->ridArray[i+mid];

        // erase the moved entries from left node
        leftNode->keyArray[i+mid] = 0;
        leftNode->ridArray[i+mid].page_number = 0;
    }

    // insert and copy up
    if (keyOps < rightNode->keyArray[0]) {
        insertLeaf(leftPage, key, rid);
    }
    else insertLeaf(rightPage, key, rid);

    // update the sibling pointers
    rightNode->rightSibPageNo = leftNode->rightSibPageNo;
    leftNode->rightSibPageNo = rightId;

    // update recurPair for the parent level
    recurPair = new PageKeyPair<int>();
    recurPair->set(rightId, rightNode->keyArray[0]);
}

// -----------------------------------------------------------------------------
// BTreeIndex::splitNonLeaf
// -----------------------------------------------------------------------------

const void BTreeIndex::splitNonLeaf(Page* curPage, PageId leafPageId, PageKeyPair<int> *&recurPair)
{
    PageId rightId;
    Page* rightPage;
    this->bufMgr->allocPage(this->file, rightId, rightPage);
    this->updateMidNode(curPage, rightPage, rightId, this->recurPair);

    // if the curNode is the root
    if (leafPageId == this->rootPageNum)
    {
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

    this->bufMgr->unPinPage(this->file, leafPageId, true);
    this->bufMgr->unPinPage(this->file, rightId, true);
}

// -----------------------------------------------------------------------------
// BTreeIndex::updateMidNode
// -----------------------------------------------------------------------------

const void BTreeIndex::updateMidNode(Page* leftPage, Page* rightPage, PageId rightId,
                  PageKeyPair<int> *&recurPair)
{
    // std::cout << "------------updateMidLeaf used------------" << std::endl;
    NonLeafNodeInt* leftNode = (NonLeafNodeInt*) leftPage;
    NonLeafNodeInt* rightNode = (NonLeafNodeInt*) rightPage;

    // determine the index of the middle key, increment by one because push up
    int mid = leafOccupancy / 2 + 1;

    // move the latter half of the keys and rids to the new leaf
    for(int i = 0; i < nodeOccupancy-mid; i++) {
        rightNode->keyArray[i] = leftNode->keyArray[i+mid];
        rightNode->pageNoArray[i] = leftNode->pageNoArray[i+mid+1];

        // erase the moved entries from left node
        leftNode->keyArray[i+mid] = 0;
        leftNode->pageNoArray[i+mid+1] = 0;
    }

    // update recurPair for the parent level
    recurPair = new PageKeyPair<int>();
    recurPair->set(rightId, leftNode->keyArray[mid-1]);

    // push up instead of copy up
    leftNode->keyArray[mid-1] = 0;
    leftNode->pageNoArray[mid] = 0;

    // insert the new child entry
    if (recurPair->key < rightNode->keyArray[0]) {
        insertNonLeaf(leftPage, recurPair);
    } 
    else insertNonLeaf(rightPage, recurPair);

    // update the level attribute
    rightNode->level = leftNode->level;
}


// -----------------------------------------------------------------------------
// BTreeIndex::findPageNoInNonLeaf
// -----------------------------------------------------------------------------

const void BTreeIndex::findPageNoInNonLeaf(NonLeafNodeInt *curNode, PageId &nextNodeNum, const void *key)
{
    if (this->attributeType==INTEGER) {
        int keyOp = *(int*) key;
        for (int i = nodeOccupancy - 1; i >= 0; i--) {
            if (curNode->pageNoArray[i] != 0) {
                if (i == 0 || curNode->keyArray[i-1] < keyOp) {
                    nextNodeNum = curNode->pageNoArray[i];
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
    if (scanExecuting) endScan();
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
        if(initialRootPageNum != rootPageNum) {
            this->findStartLeaf(lowValParm);
        }
        // find the smallest one that meet the predicate on leaf
        bool keyFound = false;
        this->findKeyLeaf(keyFound);
        // std::cout <<"--------"<<found<<"--------"<<std::endl;
        if (!keyFound) {
            bufMgr->unPinPage(file, currentPageNum, false);
            throw NoSuchKeyFoundException();
        }
    }
}

// -----------------------------------------------------------------------------
// BTreeIndex::findStartLeaf
// -----------------------------------------------------------------------------

const void BTreeIndex::findStartLeaf(const void* lowValParm) {
    if (this->attributeType==INTEGER) {
        NonLeafNodeInt* currentNode = (NonLeafNodeInt *) currentPageData;
        PageId nextPageNum;

        // reach toward the parent of a leaf
        while(currentNode->level != 1) {
            // find the child page given some predicate
            findPageNoInNonLeaf(currentNode, nextPageNum, lowValParm);
            bufMgr->unPinPage(file, currentPageNum, false);

            // proceed to the child page
            currentPageNum = nextPageNum;
            bufMgr->readPage(file, currentPageNum, currentPageData);
            currentNode = (NonLeafNodeInt*) currentPageData;
        }

        // parent of a leaf
        this->findPageNoInNonLeaf(currentNode, nextPageNum, lowValParm);
        bufMgr->unPinPage(file, currentPageNum, false);
        currentPageNum = nextPageNum;
        this->bufMgr->readPage(file, currentPageNum, currentPageData);
    }
}

// -----------------------------------------------------------------------------
// BTreeIndex::foundKey
// -----------------------------------------------------------------------------

const void BTreeIndex::findKeyLeaf(bool& keyFound) {
    // find the smallest one that meet the predicate
    LeafNodeInt* currentNode = (LeafNodeInt *) currentPageData;

    while ((currentNode->ridArray[0].page_number!=0)) {
        // search key on one page
        for(int i = 0; i < leafOccupancy; i++) {
            int key = currentNode->keyArray[i];

            // if key has been inserted
            if(currentNode->ridArray[i].page_number == 0) {
                break;
            }

            // find key
            if (((lowOp==GT && key>lowValInt) || (lowOp==GTE && key>=lowValInt)) && 
						((highOp==LT && key<highValInt) || (highOp==LTE && key<=highValInt))) {
                nextEntry = i;
                scanExecuting = true;
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
        if(currentNode->rightSibPageNo == 0) {
            return;
        }

        // move to the next leaf
        bufMgr->unPinPage(file, currentPageNum, false);
        currentPageNum = currentNode->rightSibPageNo;
        bufMgr->readPage(file, currentPageNum, currentPageData);
        currentNode = (LeafNodeInt *) currentPageData;
    }
}

// -----------------------------------------------------------------------------
// BTreeIndex::scanNext
// -----------------------------------------------------------------------------

const void BTreeIndex::scanNext(RecordId& outRid) 
{
  if(!scanExecuting)
  {
    throw ScanNotInitializedException();
  }
	// Cast page to node
	LeafNodeInt* currentNode = (LeafNodeInt *) currentPageData;
  if(currentNode->ridArray[nextEntry].page_number == 0 or nextEntry == leafOccupancy)
  {
    // Unpin page and read papge
    bufMgr->unPinPage(file, currentPageNum, false);
    // No more next leaf
    if(currentNode->rightSibPageNo == 0)
    {
      throw IndexScanCompletedException();
    }
    currentPageNum = currentNode->rightSibPageNo;
    bufMgr->readPage(file, currentPageNum, currentPageData);
    currentNode = (LeafNodeInt *) currentPageData;
    // Reset nextEntry
    nextEntry = 0;
  }
 
  // Check  if rid satisfy
  int key = currentNode->keyArray[nextEntry];
  if (((lowOp==GT && key>lowValInt) || (lowOp==GTE && key>=lowValInt)) && 
			((highOp==LT && key<highValInt) || (highOp==LTE && key<=highValInt))) {
    outRid = currentNode->ridArray[nextEntry];
    // Incrment nextEntry
    nextEntry++;
    // If current page has been scanned to its entirety
  }
  else
  {
    throw IndexScanCompletedException();
  }
}

// -----------------------------------------------------------------------------
// BTreeIndex::endScan
// -----------------------------------------------------------------------------

const void BTreeIndex::endScan() 
{
  if(!scanExecuting)
  {
    throw ScanNotInitializedException();
  }
  scanExecuting = false;
  // Unpin page
  bufMgr->unPinPage(file, currentPageNum, false);
  // Reset variable
  currentPageData = nullptr;
  currentPageNum = static_cast<PageId>(-1);
  nextEntry = -1;
}

}