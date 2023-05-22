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
		   file, headerPageNum

		2: rootPageNum
	*/

	/* VARIABLES SCANNING
		1: scanExecuting, 

		2: nextEntry, currentPageNum, currentPageData,
		lowValInt (double/string), highValInt (double/string),
		lowOp, highOp
	*/
	
	// the constructor first checks if the specified index file exists
	std::ostringstream idxStr;
	idxStr << relationName << '.' << attrByteOffset;
	std::string indexName = idxStr.str(); // indexName is the name of the index file
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
	this->headerPageNum = 1; 
	this->InitLeaf = true;
	this->onlyRoot = true;


	// initialize pointers for first page of blobfile
	BlobFile* check_index;
	IndexMetaInfo* blob_meta;

	try {
		// file exists
		check_index = new BlobFile(relationName, false);
		this->scanExecuting = false; // no need to scan
		bool flush = false;

		this->file = (File*) check_index;

		// check metapage
		Page* blob_metaPage;
		try {
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
		catch (BadIndexInfoException& e) {
			throw BadIndexInfoException("Meta info does not match");
		}
	}
	catch (FileNotFoundException& e) {
		// file does not exist -> create
		check_index = new BlobFile(relationName, true);
		this->file = (File*) check_index;

		// define the meta page info from new page in memory
		// PageNo=1
		Page* blob_metaPage;
		this->bufMgr->allocPage(this->file, this->headerPageNum, blob_metaPage);
		blob_meta = (IndexMetaInfo*) blob_metaPage;
		strncpy(blob_meta->relationName, meta_name, sizeof(blob_meta->relationName));
		blob_meta->attrByteOffset = attrByteOffset;
		blob_meta->attrType = attrType; 

		// assign rootpage id through page allocation
		this->bufMgr->allocPage(this->file, this->rootPageNum, this->rootP);
		blob_meta->rootPageNo = this->rootPageNum;

		// finish with the meta data, unpin the allocated page
		this->bufMgr->unPinPage(this->file, this->headerPageNum, true);

		if (this->attributeType==INTEGER) {
			LeafNodeInt* rootN = (LeafNodeInt*) this->rootP;
			rootN->rightSibPageNo = 0;
			rootN->slot_occupied = 0;
			for(int i = 0; i < (this->leafOccupancy); i++) {
				rootN->keyArray[i] = 0;
			}
		}
		else if (this->attributeType==DOUBLE) {
			LeafNodeDouble* rootN = (LeafNodeDouble*) this->rootP;
			rootN->rightSibPageNo = 0;
			rootN->slot_occupied = 0;
			for(int i = 0; i < (this->leafOccupancy); i++) {
				rootN->keyArray[i] = 0;
			}
		}
		else {
			LeafNodeString* rootN = (LeafNodeString*) this->rootP;
			rootN->rightSibPageNo = 0;
			rootN->slot_occupied = 0;
			for(int i = 0; i < (this->leafOccupancy); i++) {
				strncpy(rootN->keyArray[i], "", sizeof(rootN->keyArray[i]));
			}
		}

		// unpin root page
		this->bufMgr->unPinPage(this->file, this->rootPageNum, true);

		// scan file and insert all tuples 
		FileScan* fileScan = new FileScan(relationName, this->bufMgr);
		RecordId scanRid;
		while (1) {
			try {
				fileScan->scanNext(scanRid);
				std::string curRecord = fileScan->getRecord();
				if ((this->attributeType==INTEGER) || 
					(this->attributeType==DOUBLE)) {
					this->insertEntry((void*) (curRecord.c_str() + attrByteOffset), 
										scanRid);
				}
				else {
					char key[STRINGSIZE + 1];
                    char* src = (char*)(curRecord.c_str() + attrByteOffset);
                    strncpy(key, src, sizeof(key));
					std::string strKey = std::string(key);
					this->insertEntry((void*) &strKey, scanRid);
				}
			}
			catch(IndexScanCompletedException e)
			{
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

}

// -----------------------------------------------------------------------------
// BTreeIndex::insertEntry
// -----------------------------------------------------------------------------

const void BTreeIndex::insertEntry(const void *key, const RecordId rid) 
{
	if (this->attributeType==INTEGER) {
		const void *newKey;
		RecordId newRid;
		PageId lastRecurId = NULL;
		this->insertRecursive(this->rootP, this->rootPageNum, key, rid, lastRecurId);
	}
}

// -----------------------------------------------------------------------------
// BTreeIndex::insertRecursive
// -----------------------------------------------------------------------------

const void BTreeIndex::insertRecursive(Page* curPage, PageId curPageId, const void *key, 
						const RecordId rid, PageId& lastRecurId)
{
	// case when there is only one node in the tree
	if (onlyRoot) this->InitLeaf = true;

	// find root page first
	this->bufMgr->readPage(this->file, this->rootPageNum, this->rootP);

	if (this->InitLeaf) { // curPage is leaf
		this->InitLeaf = false;
		if (this->attributeType==INTEGER) {
			LeafNodeInt* leafInt = (LeafNodeInt*)curPage;
			int keyInt = *((const int*) key);

			if (leafInt->slot_occupied == leafOccupancy) { // page is full
				Page* newPage;
				PageId newPageId;
				this->splitLeaf(newPage, newPageId, curPage, 
								curPageId, key, rid, lastRecurId);
			}
			else { // page not full
				lastRecurId = NULL;
				this->insertLeaf(curPage, key, rid);
			}
		}
	}
	else { // the curPage is not a leaf
		if (this->attributeType==INTEGER) {
			NonLeafNodeInt* nonleafInt = (NonLeafNodeInt*)curPage;
			PageId nextPageId = this->findPageNoInNonLeaf(key, curPage);
			Page* nextPage;
			this->bufMgr->readPage(this->file, nextPageId, nextPage);

			// check if this node is the parent of leaf
			if (nonleafInt->level == 1) this->InitLeaf = true;

			// call on recursion
			this->insertRecursive(nextPage, nextPageId, key, rid, lastRecurId);

			// check the status of child split when the tree goes from leaf to root
			if (lastRecurId) {

				// check if need to split on this level
				if (nonleafInt->slot_occupied == leafOccupancy) { // page is full
					Page* newPage;
					PageId newPageId;
					this->splitNonLeaf(newPage, newPageId,
									   curPage, curPageId, key, lastRecurId);
				}
				else { // page is not full
					this->insertNonLeaf(curPage, key, lastRecurId);
					lastRecurId = NULL;
				}
			}
		}
	}
}

// -----------------------------------------------------------------------------
// BTreeIndex::insertLeaf
// -----------------------------------------------------------------------------

const void BTreeIndex::insertLeaf(Page* curPage, const void *key, const RecordId rid) {

	if (this->attributeType==INTEGER) {
		LeafNodeInt* leafInt = (LeafNodeInt*) curPage;
		int keyInt = *((const int*) key);
		if (leafInt->slot_occupied == 0) {
			leafInt->keyArray[0] = keyInt;
    		leafInt->ridArray[0] = rid;
			return;
		}
		int insertIdx;
		for (int i=0; i<leafInt->slot_occupied; i++) {
			if (keyInt <= leafInt->keyArray[i]) {
				insertIdx = i;
				break;
			}
			insertIdx = leafInt->slot_occupied;
		}
		
		// locate the index for the inserting key
		for (int i=(leafInt->slot_occupied); i>insertIdx; i--) {
			leafInt->keyArray[i] = leafInt->keyArray[i-1];
			leafInt->ridArray[i] = leafInt->ridArray[i-1]; 
		}

		// insert into the leaf
		leafInt->keyArray[insertIdx] = keyInt;
		leafInt->ridArray[insertIdx] = rid;
		leafInt->slot_occupied ++; // increment number of slot occupied
	}
}

// -----------------------------------------------------------------------------
// BTreeIndex::insertNonLeaf
// -----------------------------------------------------------------------------

const void BTreeIndex::insertNonLeaf(Page* curPage, const void *key, PageId pid)
{
	if (this->attributeType==INTEGER) {
		NonLeafNodeInt* nonLeafInt = (NonLeafNodeInt*) curPage;
		int keyInt = *((const int*) key);
		if (keyInt < nonLeafInt->keyArray[0]) {
			nonLeafInt->keyArray[1] = nonLeafInt->keyArray[0];
    		nonLeafInt->pageNoArray[2] = nonLeafInt->pageNoArray[1];
			nonLeafInt->keyArray[0] = keyInt;
			nonLeafInt->pageNoArray[1] = pid;
			return;
		}
		int insertIdx;
		for (int i=nonLeafInt->slot_occupied; i>0; i--) {
			if (keyInt >= nonLeafInt->keyArray[i]) {
				insertIdx = i;
				break;
			}
			nonLeafInt->keyArray[i] = nonLeafInt->keyArray[i-1];
			nonLeafInt->pageNoArray[i+1] = nonLeafInt->pageNoArray[i];
			insertIdx = i;
		}

		// insert into the leaf
		nonLeafInt->keyArray[insertIdx] = keyInt;
		nonLeafInt->pageNoArray[insertIdx+1] = pid;
		nonLeafInt->slot_occupied ++; // increment number of slot occupied
	}
}

// -----------------------------------------------------------------------------
// BTreeIndex::splitLeaf
// -----------------------------------------------------------------------------

const void BTreeIndex::splitLeaf(Page* newPage, PageId newPageId, Page* curPage, 
						PageId curPageId, const void *key, const RecordId rid, 
						PageId& lastRecurId) 
{
	this->onlyRoot = false;
	// allocate new page to a new node
	this->bufMgr->allocPage(this->file, newPageId, newPage);

	if (this->attributeType==INTEGER) {
		LeafNodeInt* origLeaf = (LeafNodeInt*) curPage;
		LeafNodeInt* newLeaf = (LeafNodeInt*) newPage;
		int keyInt = *((const int*) key);
		int splitFactor = (origLeaf->slot_occupied) / 2;
		if ((leafOccupancy % 2) == 1 && keyInt > origLeaf->keyArray[splitFactor])
		{
			splitFactor ++;
		}

		// connect the sibling leaf nodes
		newLeaf->rightSibPageNo = origLeaf->rightSibPageNo;
    	origLeaf->rightSibPageNo = newPageId;

		// redistribute between two nodes
		int newSlot, origSlot;
		origSlot = origLeaf->slot_occupied;
		newSlot = 0;
		for (int i=0; i<origLeaf->slot_occupied-splitFactor; i++) {
			newLeaf->keyArray[i] = origLeaf->keyArray[i+splitFactor];
			newLeaf->ridArray[i] = origLeaf->ridArray[i+splitFactor];
			newSlot ++;
			origSlot --;
		}

		// recalculate the slot occupied variable
		origLeaf->slot_occupied = origSlot;
		newLeaf->slot_occupied = newSlot;

		// insert new key into either node
		if (keyInt < newLeaf->keyArray[0]) {
			this->insertLeaf(curPage, key, rid);
		}
		else {
			this->insertLeaf(newPage, key, rid);
		}

		// update the split status
		lastRecurId = newPageId;
		this->InitLeaf = false;

		// make a new root node, default is the smallest right child
		// this will happen only if the current split is at the root
		if (curPageId==this->rootPageNum) {
			Page* newRoot;
			PageId newRootId;
			this->bufMgr->allocPage(this->file, newRootId, newRoot);

			// update the new root
			NonLeafNodeInt* newRootN = (NonLeafNodeInt*) newRoot;
			newRootN->level = 1; // because this is originally leaf
			newRootN->pageNoArray[0] = curPageId;
			newRootN->pageNoArray[1] = newPageId;
			newRootN->keyArray[0] = newLeaf->keyArray[0];
			newRootN->slot_occupied = 1;
			this->rootPageNum = newRootId;

			// update meta page
			Page* headerP;
			this->bufMgr->readPage(this->file, this->headerPageNum, headerP);
			IndexMetaInfo* metaP = (IndexMetaInfo*) headerP;
			metaP->rootPageNo = newRootId;

			// free buffer memory
			this->bufMgr->unPinPage(this->file, this->headerPageNum, true);
			this->bufMgr->unPinPage(this->file, newRootId, true);
		}

		// free buffer memory
		this->bufMgr->unPinPage(this->file, curPageId, true);
  		this->bufMgr->unPinPage(this->file, newPageId, true);
	}
}

// -----------------------------------------------------------------------------
// BTreeIndex::splitNonLeaf
// -----------------------------------------------------------------------------

const void BTreeIndex::splitNonLeaf(Page* newPage, PageId newPageId, Page* curPage, 
						PageId curPageId, const void *key, PageId& lastRecurId)
{
	this->bufMgr->allocPage(this->file, newPageId, newPage);
	if (this->attributeType==INTEGER) {
		NonLeafNodeInt* origNodeInt = (NonLeafNodeInt*) curPage;
		NonLeafNodeInt* newNodeInt = (NonLeafNodeInt*)newPage;
		newNodeInt->level = origNodeInt->level;
		int keyInt = *((const int*) key);
		int splitFactor = (origNodeInt->slot_occupied) / 2;
		int pushupIndex = splitFactor;

		if (nodeOccupancy % 2 == 0)
		{
			pushupIndex = keyInt < origNodeInt->keyArray[splitFactor] ? splitFactor -1 : splitFactor;
		}
		
		splitFactor = pushupIndex + 1;
		// redistribute between two nodes
		int newSlot, origSlot;
		origSlot = origNodeInt->slot_occupied;
		newSlot = 0;
		for (int i=0; i<origNodeInt->slot_occupied-splitFactor; i++) {
			newNodeInt->keyArray[i] = origNodeInt->keyArray[i+splitFactor];
			newNodeInt->pageNoArray[i] = origNodeInt->pageNoArray[i+splitFactor+1];
			newSlot ++;
			origSlot --;
		}

		// recalculate the slot occupied variable
		origNodeInt->slot_occupied = origSlot-1;
		newNodeInt->slot_occupied = newSlot;

		// remove the entry that is pushed up from current node
		origNodeInt->keyArray[pushupIndex] = 0;
		origNodeInt->pageNoArray[pushupIndex] = 0;

		// insert new key into either node
		if (keyInt < newNodeInt->keyArray[0]) {
			this->insertNonLeaf(curPage, key, lastRecurId);
		}
		else {
			this->insertNonLeaf(newPage, key, lastRecurId);
		}

		// update the split status
		lastRecurId = newPageId;
		this->InitLeaf = false;

		// make a new root node, default is the smallest right child
		// this will happen only if the current split is at the root
		if (curPageId==this->rootPageNum) {
			Page* newRoot;
			PageId newRootId;
			this->bufMgr->allocPage(this->file, newRootId, newRoot);

			// update the new root
			NonLeafNodeInt* newRootN = (NonLeafNodeInt*) newRoot;
			newRootN->level = 999; // because this is a nonleaf to split
			newRootN->pageNoArray[0] = curPageId;
			newRootN->pageNoArray[1] = newPageId;
			newRootN->keyArray[0] = newNodeInt->keyArray[0];
			newRootN->slot_occupied = 1;
			this->rootPageNum = newRootId;

			// update meta page
			Page* headerP;
			this->bufMgr->readPage(this->file, this->headerPageNum, headerP);
			IndexMetaInfo* metaP = (IndexMetaInfo*) headerP;
			metaP->rootPageNo = newRootId;

			// free buffer memory
			this->bufMgr->unPinPage(this->file, this->headerPageNum, true);
			this->bufMgr->unPinPage(this->file, newRootId, true);
		}

		// free buffer memory
		this->bufMgr->unPinPage(this->file, curPageId, true);
  		this->bufMgr->unPinPage(this->file, newPageId, true);
	}
} 

// -----------------------------------------------------------------------------
// BTreeIndex::findPageNoInNonLeaf
// -----------------------------------------------------------------------------

PageId BTreeIndex::findPageNoInNonLeaf(const void *key, Page* page) // check its validity
{
	if (this->attributeType==INTEGER) {
		int keyInt = *((const int*) key);
		NonLeafNodeInt* node = (NonLeafNodeInt*) page;
		for (int i=0; i<node->slot_occupied; i++) {
			if (keyInt<node->keyArray[i]) {
				return node->pageNoArray[i];
			}
		}
		return node->pageNoArray[node->slot_occupied];
	}
	else if (this->attributeType==DOUBLE) {
		double keyDouble = *((const double*) key);
		NonLeafNodeDouble* node = (NonLeafNodeDouble*) page;
		for (int i=0; i<node->slot_occupied; i++) {
			if (keyDouble<node->keyArray[i]) {
				return node->pageNoArray[i];
			}
		}
		return node->pageNoArray[node->slot_occupied];
	}
	else {
		std::string keyStr = *((std::string*) key);
		NonLeafNodeString* node = (NonLeafNodeString*) page;
		for (int i=0; i<node->slot_occupied; i++) {
			if (keyStr.compare(node->keyArray[i])<0) {
				return node->pageNoArray[i];
			}
		}
		return node->pageNoArray[node->slot_occupied];
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
	if ((lowOpParm != GT || lowOpParm != GTE)
		|| (highOpParm != LT || highOpParm != LTE)) {
		throw BadOpcodesException();
	}
	if (lowValInt > highValInt) throw BadScanrangeException();
	this->lowOp = lowOpParm;
	this->highOp = highOpParm;
	this->currentPageNum = this->rootPageNum;
	this->bufMgr->readPage(this->file, this->currentPageNum, this->currentPageData);

	if (this->attributeType==INTEGER) {
		this->lowValInt = *(int*) lowValParm;
        this->highValInt = *(int*) highValParm;

		if(!(this->onlyRoot)) { // If the root is not a leaf
			NonLeafNodeInt* curNode = (NonLeafNodeInt*) this->currentPageData;
			PageId nextPageNum;

			// Loop until we find a leaf
			while(curNode->level != 1)
			{
				PageId nextId = this->findPageNoInNonLeaf(lowValParm, this->currentPageData);
				this->bufMgr->unPinPage(this->file, this->currentPageNum, false); 
				this->currentPageNum = nextId; 
				this->bufMgr->readPage(this->file, this->currentPageNum, this->currentPageData); 
				curNode = (NonLeafNodeInt*)(this->currentPageData);
			}

			// now curNode is the parent of leaf node
			PageId nextId = this->findPageNoInNonLeaf(lowValParm, this->currentPageData);
			this->bufMgr->unPinPage(this->file, this->currentPageNum, false);
			this->bufMgr->readPage(this->file, this->currentPageNum, this->currentPageData);
		}
		LeafNodeInt* curLeaf = (LeafNodeInt*) this->currentPageData; // leaf node
		
		// check if the page is empty
		if (curLeaf->slot_occupied == 0){
			throw NoSuchKeyFoundException();
		}

		// find the key
		while (1) {
			// Search from the left leaf page to the right to find the fit
			bool nullVal = false;
			for(int i = 0; i < leafOccupancy && !nullVal; i++)
			{
			int key = curLeaf->keyArray[i];
			// Check if the next one in the key is not inserted
			if(i < leafOccupancy - 1 && i > curLeaf->slot_occupied)
			{
				nullVal = true;
			}
			
			if(checkKey(lowValInt, lowOp, highValInt, highOp, key))
			{
				// select
				this->nextEntry = i;
				this->scanExecuting = true;
				break;
			}
			else if((this->highOp == LT && key >= highValInt) || 
					(this->highOp == LTE && key > highValInt))
			{
				bufMgr->unPinPage(this->file, this->currentPageNum, false);
				throw NoSuchKeyFoundException();
			}
			
			// Did not find any matching key in this leaf, go to next leaf
			if(i == leafOccupancy - 1 || nullVal){
				//unpin page
				this->bufMgr->unPinPage(this->file, this->currentPageNum, false);
				//did not find the matching one in the more right leaf
				if(curLeaf->rightSibPageNo == 0) {
					throw NoSuchKeyFoundException();
				}
				this->currentPageNum = curLeaf->rightSibPageNo;
				this->bufMgr->readPage(this->file, this->currentPageNum, currentPageData);
			}
			}
		}
	}
	else if (this->attributeType==DOUBLE) {
		this->lowValInt = *(double*)lowValParm;
        this->highValInt = *(double*) highValParm;
	}
	else {
		this->lowValString = *(std::string*)lowValParm;
		this->highValString = *(std::string*)highValParm;
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
	if (this->attributeType==INTEGER) {
		LeafNodeInt* currentNode = (LeafNodeInt *) this->currentPageData;
		if(nextEntry==currentNode->slot_occupied || nextEntry == leafOccupancy)
		{
			// Unpin page and read papge
			this->bufMgr->unPinPage(this->file, this->currentPageNum, false);
			// No more next leaf
			if(currentNode->rightSibPageNo == 0)
			{
				throw IndexScanCompletedException();
			}
			currentPageNum = currentNode->rightSibPageNo;
			this->bufMgr->readPage(this->file, this->currentPageNum, this->currentPageData);
			currentNode = (LeafNodeInt *) this->currentPageData;
			// Reset nextEntry
			this->nextEntry = 0;
		}

		// Check  if rid satisfy
		int key = currentNode->keyArray[nextEntry];
		if(checkKey(lowValInt, lowOp, highValInt, highOp, key))
		{
			outRid = currentNode->ridArray[nextEntry];
			// Incrment nextEntry
			nextEntry++;
		}
		else
		{
			throw IndexScanCompletedException();
		}
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
	this->scanExecuting = false;
	// Unpin page
	bufMgr->unPinPage(this->file, this->currentPageNum, false);
	// Reset variable
	this->currentPageData = nullptr;
	this->currentPageNum = PageId (-1);
	this->nextEntry = -1;
}

// -----------------------------------------------------------------------------
// BTreeIndex::checkKey
// -----------------------------------------------------------------------------

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
