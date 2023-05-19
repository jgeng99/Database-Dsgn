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
	this->headerPageNum = 1; 
	this->rootIsLeaf = true;


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
		PageId blob_metaPageId;
		this->bufMgr->allocPage(this->file, blob_metaPageId, blob_metaPage);
		blob_meta = (IndexMetaInfo*) blob_metaPage;
		strncpy(blob_meta->relationName, meta_name, sizeof(blob_meta->relationName));
		blob_meta->attrByteOffset = attrByteOffset;
		blob_meta->attrType = attrType; 

		// assign rootpage id through page allocation
		this->bufMgr->allocPage(this->file, this->rootPageNum, this->rootP);
		blob_meta->rootPageNo = this->rootPageNum;

		// finish with the meta data, unpin the allocated page
		this->bufMgr->unPinPage(this->file, this->headerPageNum, true);

		// update rootpagenum
		this->rootPageNum = this->headerPageNum + 1;

		// create root of the tree because file does not exist
		this->bufMgr->allocPage(this->file, this->rootPageNum, this->rootP);
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
			for(int i = 0; i < (this->leafOccupancy); i++) {
				strncpy(rootN->keyArray[i], "", sizeof(rootN->keyArray[i]));
			}
		}

		// scan file and insert all tuples 
		FileScan* fileScan = new FileScan(relationName, this->bufMgr);
		this->scanExecuting = true;
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
	if (this->rootIsLeaf) {
		this->rootIsLeaf = false;
		if (this->attributeType==INTEGER) {
			this->insertLeaf(key, rid);
		}
	}
}

// -----------------------------------------------------------------------------
// BTreeIndex::insertLeaf
// -----------------------------------------------------------------------------

const void BTreeIndex::insertLeaf(const void *key, const RecordId rid) {
	// find root page first
	this->bufMgr->readPage(this->file, this->rootPageNum, this->rootP);

	if (this->attributeType==INTEGER) {
		LeafNodeInt* leafInt = (LeafNodeInt*) this->rootP;
		int keyInt = *((const int*) key);
		if (leafInt->slot_occupied == leafOccupancy) { // page is full
			Page* newPage;
			PageId newPageId;
			this->splitLeaf(newPage, newPageId, this->rootP, key, rid);
            makeNewRootNode(rootPageNum, newPage, true);
		}
		else { // page not full
			int insertIdx;
			for (int i=0; i<leafInt->slot_occupied; i++) {
				if (keyInt <= leafInt->keyArray[i]) {
					insertIdx = i;
					break;
				}
			}
			
			// locate the index for the inserting key
			for (int i=(leafInt->slot_occupied-1); i>insertIdx; i--) {
				leafInt->keyArray[i] = leafInt->keyArray[i-1];
				leafInt->ridArray[i] = leafInt->ridArray[i-1];
			}

			// insert into the leaf
			leafInt->keyArray[insertIdx] = keyInt;
			leafInt->ridArray[insertIdx] = rid;
			leafInt->slot_occupied ++; // increment number of slot occupied
		}
	}
}

// -----------------------------------------------------------------------------
// BTreeIndex::splitLeaf
// -----------------------------------------------------------------------------

const void BTreeIndex::splitLeaf(Page* newPage, PageId newPageId, Page* curPage, 
						const void *key, const RecordId rid) 
{
	// allocate new page to a new node
	this->bufMgr->allocPage(this->file, newPageId, newPage);
	if (this->attributeType==INTEGER) {
		LeafNodeInt* origLeaf = (LeafNodeInt*) curPage;
		LeafNodeInt* newLeaf = (LeafNodeInt*) newPage;
		int splitFactor = (origLeaf->slot_occupied) / 2;

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

		// 
	}
}

// -----------------------------------------------------------------------------
// BTreeIndex::insertRecursive
// -----------------------------------------------------------------------------

const void BTreeIndex::insertRecursive(const void *key, const RecordId rid, 
									Page* page) 
{
	if (this->attributeType==INTEGER) {
		NonLeafNodeInt* node = (NonLeafNodeInt*) page;

		// check if this is the first time being inserted
		if (node->keyArray[0]==0) {
			node->keyArray[0] = *((const int*) key);
			node->level ++;

			// create two leaf node, one left and one right, to the current root
			Page* lleafP;
			PageId lleafNum;
			bufMgr->allocPage(this->file, lleafNum, lleafP);
			LeafNodeInt* lleafN = (LeafNodeInt*) lleafP;

			Page* rleafP; 
			PageId rleafNum;
			bufMgr->allocPage(file, rleafNum, rleafP);
			LeafNodeInt* rleafN = (LeafNodeInt*) rleafP;
			
			// link root and leaf page
			node->pageNoArray[0] = lleafNum;
			node->pageNoArray[1] = rleafNum;

			// fill the array with initial values
			for(int i = 0; i < leafOccupancy; i++) {
				lleafN->keyArray[i] = 0;
				rleafN->keyArray[i] = 0;
			}

			// link the leaf nodes
			lleafN->rightSibPageNo = rleafNum;
			rleafN->rightSibPageNo = NULL;

			// set the right leaf node to have the inserted value
			rleafN->keyArray[0] = *((const int*) key);
			rleafN->ridArray[0] = rid;

			bufMgr->unPinPage(this->file, lleafNum, true);
			bufMgr->unPinPage(this->file, rleafNum, true);

			// if nonleaf node
			if (...) {
				PageId curId = findPageNoInNonLeaf(key, page);
				Page *childP;
				this->bufMgr->readPage(this->file, curId, childP);
				insertRecursive(key, rid, childP);
			}
			// if leaf node
			else {

			}
		}
	}
}

// -----------------------------------------------------------------------------
// BTreeIndex::findPageNoInNonLeaf
// -----------------------------------------------------------------------------

PageId BTreeIndex::findPageNoInNonLeaf(const void *key, Page* page) 
{
	if (this->attributeType==INTEGER) {
		int int_key = *((const int*) key);
		NonLeafNodeInt* node = (NonLeafNodeInt*) page;
		int i = 0;
		while (i < leafOccupancy) {
			i ++;
			if (i == leafOccupancy) {
				break;
			}
			if (node->keyArray[0] > int_key) {
				return node->pageNoArray[0];
			}
			else if ((int_key >= node->keyArray[i-1]) 
					&& (int_key < node->keyArray[i])) {
				return node->pageNoArray[i];
			}
		}
		
		return node->pageNoArray[leafOccupancy];
	}
	else if (this->attributeType==DOUBLE) {
		double double_key = *((const double*) key);
		NonLeafNodeDouble* node = (NonLeafNodeDouble*) page;
		int i = 0;
		while (i < leafOccupancy) {
			i ++;
			if (i == leafOccupancy) {
				break;
			}
			if (node->keyArray[0] > double_key) {
				return node->pageNoArray[0];
			}
			else if ((double_key >= node->keyArray[i-1]) 
					&& (double_key < node->keyArray[i])) {
				return node->pageNoArray[i];
			}
		}
		
		return node->pageNoArray[leafOccupancy];
	}
	else {
		std::string str_key = *((std::string*) key);
		NonLeafNodeString* node = (NonLeafNodeString*) page;
		int i = 0;
		while (i < leafOccupancy) {
			i++;
			if (i == leafOccupancy) {
				break;
			}
			if(strcmp(str_key.c_str(), node->keyArray[0]) < 0) {
				return node->pageNoArray[0]; 
			}
			else if ((strcmp(str_key.c_str(), node->keyArray[i-1]) >= 0)
					&& (strcmp(str_key.c_str(), node->keyArray[i]) < 0)) {
				return node->pageNoArray[1];
			}
		}

		return node->pageNoArray[leafOccupancy];
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

}

// -----------------------------------------------------------------------------
// BTreeIndex::scanNext
// -----------------------------------------------------------------------------

const void BTreeIndex::scanNext(RecordId& outRid) 
{

}

// -----------------------------------------------------------------------------
// BTreeIndex::endScan
// -----------------------------------------------------------------------------
//
const void BTreeIndex::endScan() 
{

}

}
