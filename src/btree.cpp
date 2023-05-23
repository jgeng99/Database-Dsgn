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
    this->InitLeaf = true;
    this->onlyRoot = true;

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
  RIDKeyPair<int> dataEntry;
  dataEntry.set(rid, *((int *)key));
  // root
  Page* root;
  // PageId rootPageNum;
  bufMgr->readPage(file, rootPageNum, root);
  recurPair = nullptr;

  insertRecursive(root, rootPageNum, initialRootPageNum == rootPageNum ? true : false, 
                  key, rid, recurPair);
}

// -----------------------------------------------------------------------------
// BTreeIndex::insertEntry
// -----------------------------------------------------------------------------

const void BTreeIndex::insertRecursive(Page *curPage, PageId curPageNum, bool nodeIsLeaf, 
                const void *key, const RecordId rid, PageKeyPair<int> *&recurPair)
{
  if (nodeIsLeaf)
  {
    // LeafNodeInt* leafInt = (LeafNodeInt*) curPage;
    // page is not full
    if (!leafFull(curPage))
    {
      insertLeaf(curPage, key, rid);
      bufMgr->unPinPage(file, curPageNum, true);
      recurPair = nullptr;
    }
    else
    {
      splitLeaf(curPage, curPageNum, key, rid, recurPair);
    }
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
    
    // no split in child, just return
    if (recurPair == nullptr)
    {
	    // unpin current page from call stack
	    bufMgr->unPinPage(file, curPageNum, false);
    }
    else
	  { 
      // if the curpage is not full
      if (!nodeFull(curPage))
      {
        // insert the newchildEntry to curpage
        insertNonLeaf(curPage, recurPair);
        recurPair = nullptr;
        // finish the insert process, unpin current page
        bufMgr->unPinPage(file, curPageNum, true);
      }
      else
      {
        splitNonLeaf(curPage, curPageNum, key, recurPair);
      }
    }
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
// BTreeIndex::leafFull
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

const void BTreeIndex::insertNonLeaf(Page* curPage, PageKeyPair<int> *entry)
{
  NonLeafNodeInt* nonLeaf = (NonLeafNodeInt*) curPage;

  // determine the number of valid keys in the non-leaf node
  int validKeys = 0;
  while((validKeys<nodeOccupancy) && 
      (nonLeaf->pageNoArray[validKeys+1]!=0)) {
    validKeys++;
  }

  // find insertion point
  int start = 0, end = validKeys - 1;
  while(start <= end) {
    // use binary search to save runtime because sorted
    int mid = start + (end - start) / 2;
    if(nonLeaf->keyArray[mid] < entry->key) {
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
  nonLeaf->keyArray[insertionPoint] = entry->key;
  nonLeaf->pageNoArray[insertionPoint + 1] = entry->pageNo;
}

// -----------------------------------------------------------------------------
// BTreeIndex::splitLeaf
// -----------------------------------------------------------------------------

const void BTreeIndex::splitLeaf(Page* curPage, PageId leafPageNum, 
          const void *key, const RecordId rid, PageKeyPair<int> *&newchildEntry)
{
  LeafNodeInt* leafInt = (LeafNodeInt*) curPage;
  int keyInt = *(int*) key;
  PageId newPageNum;
  Page *newPage;
  bufMgr->allocPage(file, newPageNum, newPage);
  LeafNodeInt *newLeafNode = (LeafNodeInt *)newPage;

  int mid = (leafOccupancy)/2;
  // odd number of keys
  if (leafOccupancy %2 == 1 && keyInt > leafInt->keyArray[mid])
  {
    mid = mid + 1;
  }
  // copy half the page to newLeafNode
  for(int i = mid; i < leafOccupancy; i++)
  {
    newLeafNode->keyArray[i-mid] = leafInt->keyArray[i];
    newLeafNode->ridArray[i-mid] = leafInt->ridArray[i];
    leafInt->keyArray[i] = 0;
    leafInt->ridArray[i].page_number = 0;
  }

  
  if (keyInt < leafInt->keyArray[mid-1])
  {
    insertLeaf(curPage, key, rid);
  }
  else
  {
    insertLeaf(newPage, key, rid);
  }

  // update sibling pointer
  newLeafNode->rightSibPageNo = leafInt->rightSibPageNo;
  leafInt->rightSibPageNo = newPageNum;

  // the smallest key from second page as the new child entry
  newchildEntry = new PageKeyPair<int>();
  PageKeyPair<int> newKeyPair;
  newKeyPair.set(newPageNum, newLeafNode->keyArray[0]);
  newchildEntry = &newKeyPair;
  bufMgr->unPinPage(file, leafPageNum, true);
  bufMgr->unPinPage(file, newPageNum, true);

  // if curr page is root
  if (leafPageNum == rootPageNum)
  {
    // create a new root 
    PageId newRootPageNum;
    Page *newRoot;
    bufMgr->allocPage(file, newRootPageNum, newRoot);
    NonLeafNodeInt *newRootPage = (NonLeafNodeInt *)newRoot;

    // update metadata
    newRootPage->level = 1;
    newRootPage->pageNoArray[0] = leafPageNum;
    newRootPage->pageNoArray[1] = newPageNum;
    newRootPage->keyArray[0] = newLeafNode->keyArray[0];

    Page *meta;
    bufMgr->readPage(file, headerPageNum, meta);
    IndexMetaInfo *metaPage = (IndexMetaInfo *)meta;
    metaPage->rootPageNo = newRootPageNum;
    this->rootPageNum = newRootPageNum;
    // unpin unused page
    bufMgr->unPinPage(file, headerPageNum, true);
    bufMgr->unPinPage(file, newRootPageNum, true);
  }
}

// -----------------------------------------------------------------------------
// BTreeIndex::splitNonLeaf
// -----------------------------------------------------------------------------

const void BTreeIndex::splitNonLeaf(Page* curPage, PageId curPageId,
                      const void *key, PageKeyPair<int> *&newchildEntry)
{
  NonLeafNodeInt* nonLeafInt = (NonLeafNodeInt*) curPage;
  int keyInt = *(int*) key;
  PageId newPageNum;
  Page *newPage;
  bufMgr->allocPage(file, newPageNum, newPage);
  NonLeafNodeInt *newNode = (NonLeafNodeInt *)newPage;

  int mid = nodeOccupancy/2;
  int pushupIndex = mid;
  PageKeyPair<int> pushupEntry;
  // even number of keys
  if (nodeOccupancy % 2 == 0)
  {
    pushupIndex = newchildEntry->key < nonLeafInt->keyArray[mid] ? mid -1 : mid;
  }
  pushupEntry.set(newPageNum, nonLeafInt->keyArray[pushupIndex]);

  mid = pushupIndex + 1;
  // move half the entries to the new node
  for(int i = mid; i < nodeOccupancy; i++)
  {
    newNode->keyArray[i-mid] = nonLeafInt->keyArray[i];
    newNode->pageNoArray[i-mid] = nonLeafInt->pageNoArray[i+1];
    nonLeafInt->pageNoArray[i+1] = (PageId) 0;
    nonLeafInt->keyArray[i+1] = 0;
  }

  newNode->level = nonLeafInt->level;
  // remove the entry that is pushed up from current node
  nonLeafInt->keyArray[pushupIndex] = 0;
  nonLeafInt->pageNoArray[pushupIndex] = (PageId) 0;
  // insert the new child entry
  if (newchildEntry->key < newNode->keyArray[0]) insertNonLeaf(curPage, newchildEntry);
  else insertNonLeaf(newPage, newchildEntry);
  newchildEntry = &pushupEntry;
  bufMgr->unPinPage(file, curPageId, true);
  bufMgr->unPinPage(file, newPageNum, true);

  // if the curNode is the root
  if (curPageId == this->rootPageNum)
  {
    // create a new root 
    PageId newRootPageNum;
    Page *newRoot;
    bufMgr->allocPage(file, newRootPageNum, newRoot);
    NonLeafNodeInt *newRootPage = (NonLeafNodeInt *)newRoot;

    // update metadata
    newRootPage->level = 0;
    // newRootPage->level = 0;
    newRootPage->pageNoArray[0] = curPageId;
    newRootPage->pageNoArray[1] = newchildEntry->pageNo;
    newRootPage->keyArray[0] = newchildEntry->key;

    Page *meta;
    bufMgr->readPage(file, headerPageNum, meta);
    IndexMetaInfo *metaPage = (IndexMetaInfo *)meta;
    metaPage->rootPageNo = newRootPageNum;
    rootPageNum = newRootPageNum;
    // unpin unused page
    bufMgr->unPinPage(file, headerPageNum, true);
    bufMgr->unPinPage(file, newRootPageNum, true);
  }
}


// -----------------------------------------------------------------------------
// BTreeIndex::findPageNoInNonLeaf
// -----------------------------------------------------------------------------

const void BTreeIndex::findPageNoInNonLeaf(NonLeafNodeInt *curNode, PageId &nextNodeNum, const void *key)
{
    int keyIny = *(int*) key;
    for (int i = nodeOccupancy - 1; i >= 0; i--)
    {
        if (curNode->pageNoArray[i] != 0)
        {
            if (i == 0 || curNode->keyArray[i-1] < keyIny)
            {
                nextNodeNum = curNode->pageNoArray[i];
                return;
            }
        }
    }
}

/**
 * Begin a filtered scan of the index.  For instance, if the method is called 
 * using ("a",GT,"d",LTE) then we should seek all entries with a value 
 * greater than "a" and less than or equal to "d".
 * If another scan is already executing, that needs to be ended here.
 * Set up all the variables for scan. Start from root to find out the leaf page that contains the first RecordID
 * that satisfies the scan parameters. Keep that page pinned in the buffer pool.
 * @param lowVal  Low value of range, pointer to integer / double / char string
 * @param lowOp   Low operator (GT/GTE)
 * @param highVal High value of range, pointer to integer / double / char string
 * @param highOp  High operator (LT/LTE)
 * @throws  BadOpcodesException If lowOp and highOp do not contain one of their their expected values 
 * @throws  BadScanrangeException If lowVal > highval
 * @throws  NoSuchKeyFoundException If there is no key in the B+ tree that satisfies the scan criteria.
**/

const void BTreeIndex::startScan(const void* lowValParm,
           const Operator lowOpParm,
           const void* highValParm,
           const Operator highOpParm)
{
  
  lowValInt = *((int *)lowValParm);
  highValInt = *((int *)highValParm);

  if(!((lowOpParm == GT or lowOpParm == GTE) and (highOpParm == LT or highOpParm == LTE)))
  {
    throw BadOpcodesException();
  }
  if(lowValInt > highValInt)
  {
    throw BadScanrangeException();
  }

  lowOp = lowOpParm;
  highOp = highOpParm;

  // Scan is already started
  if(scanExecuting)
  {
    endScan();
  }

  currentPageNum = rootPageNum;
  // Start scanning by reading rootpage into the buffer pool
  bufMgr->readPage(file, currentPageNum, currentPageData);

  // root is not a leaf
  if(initialRootPageNum != rootPageNum)
  {
    // Cast
    NonLeafNodeInt* currentNode = (NonLeafNodeInt *) currentPageData;
    bool foundLeaf = false;
    while(!foundLeaf)
    {
      // Cast page to node
      currentNode = (NonLeafNodeInt *) currentPageData;
      // Check if this is the level above the leaf, if yes, the next level is the leaf
      if(currentNode->level == 1)
      {
        foundLeaf = true;
      }

      // Find the leaf
      PageId nextPageNum;
      // std::cout << "startscan: "<< currentNode->slot_occupied<<std::endl;
      findPageNoInNonLeaf(currentNode, nextPageNum, lowValParm);
      // Unpin
      bufMgr->unPinPage(file, currentPageNum, false);
      currentPageNum = nextPageNum;
      // read the nextPage
      bufMgr->readPage(file, currentPageNum, currentPageData);
    }
  }
  // Now the curNode is leaf node try to find the smallest one that satisefy the OP
  bool found = false;
  while(!found){
    // Cast page to node
    LeafNodeInt* currentNode = (LeafNodeInt *) currentPageData;
    // Check if the whole page is null
    if(currentNode->ridArray[0].page_number == 0)
    {
      bufMgr->unPinPage(file, currentPageNum, false);
      throw NoSuchKeyFoundException();
    }
    // Search from the left leaf page to the right to find the fit
    bool nullVal = false;
    for(int i = 0; i < leafOccupancy and !nullVal; i++)
    {
      int key = currentNode->keyArray[i];
      // Check if the next one in the key is not inserted
      if(i < leafOccupancy - 1 and currentNode->ridArray[i + 1].page_number == 0)
      {
        nullVal = true;
      }
      
      if(checkKey(lowValInt, lowOp, highValInt, highOp, key))
      {
        // select
        nextEntry = i;
        found = true;
        scanExecuting = true;
        break;
      }
      else if((highOp == LT and key >= highValInt) or (highOp == LTE and key > highValInt))
      {
        bufMgr->unPinPage(file, currentPageNum, false);
        throw NoSuchKeyFoundException();
      }
      
      // Did not find any matching key in this leaf, go to next leaf
      if(i == leafOccupancy - 1 or nullVal){
        //unpin page
        bufMgr->unPinPage(file, currentPageNum, false);
        //did not find the matching one in the more right leaf
        if(currentNode->rightSibPageNo == 0)
        {
          throw NoSuchKeyFoundException();
        }
        currentPageNum = currentNode->rightSibPageNo;
        bufMgr->readPage(file, currentPageNum, currentPageData);
      }
    }
  }
}

/**
  * Fetch the record id of the next index entry that matches the scan.
  * Return the next record from current page being scanned. If current page has been scanned to its entirety, move on to the right sibling of current page, if any exists, to start scanning that page. Make sure to unpin any pages that are no longer required.
  * @param outRid RecordId of next record found that satisfies the scan criteria returned in this
  * @throws ScanNotInitializedException If no scan has been initialized.
  * @throws IndexScanCompletedException If no more records, satisfying the scan criteria, are left to be scanned.
**/
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
  if(checkKey(lowValInt, lowOp, highValInt, highOp, key))
  {
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

/**
  * Terminate the current scan. Unpin any pinned pages. Reset scan specific variables.
  * @throws ScanNotInitializedException If no scan has been initialized.
**/
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

/**
  * Helper function to check if the key is satisfies
  * @param lowVal   Low value of range, pointer to integer / double / char string
  * @param lowOp    Low operator (GT/GTE)
  * @param highVal  High value of range, pointer to integer / double / char string
  * @param highOp   High operator (LT/LTE)
  * @param val      Value of the key
  * @return True if satisfies False if not
  *
**/
const bool BTreeIndex::checkKey(int lowVal, const Operator lowOp, int highVal, const Operator highOp, int key)
{
  if(lowOp == GTE && highOp == LTE)
  {
    return key <= highVal && key >= lowVal;
  }
  else if(lowOp == GT && highOp == LTE)
  {
    return key <= highVal && key > lowVal;
  }
  else if(lowOp == GTE && highOp == LT)
  {
    return key < highVal && key >= lowVal;
  }
  else
  {
    return key < highVal && key > lowVal;
  }
}

}