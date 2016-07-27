#include "utils.h"

static void fillValues(bst_t *botTree, const uint64_t botVal, const uint64_t tau){

	botTree->botVal = botVal;
	
	if (botVal > UINT64_MAX-tau)
		botTree->topVal = tau-(UINT64_MAX-botVal);
	else
		botTree->topVal = botVal + tau;

}


static bst_t* insertCenter(uint32_t newBot, uint32_t newTop, const uint64_t tau, const uint64_t *bot, bst_t *botTree){

	uint64_t newCenter;


	if (newBot == newTop){
		fillValues(&(botTree[newBot]), bot[newBot], tau);
		botTree[newBot].leftPtr = NULL;
		botTree[newBot].rightPtr = NULL;
		return &(botTree[newBot]);
	}

	if (newTop == newBot+1){
		fillValues(&(botTree[newBot]), bot[newBot], tau);
		fillValues(&(botTree[newTop]), bot[newTop], tau);
		botTree[newBot].leftPtr = NULL;
		botTree[newBot].rightPtr = &(botTree[newTop]);
		botTree[newTop].leftPtr = NULL;
		botTree[newTop].rightPtr = NULL;
		return &(botTree[newBot]);
	}

	newCenter = newBot + (newTop - newBot)/2;
	fillValues(&(botTree[newCenter]), bot[newCenter], tau);

	// left
	botTree[newCenter].leftPtr = insertCenter(newBot, newCenter-1, tau, bot, botTree);
	// right
	botTree[newCenter].rightPtr = insertCenter(newCenter+1, newTop, tau, bot, botTree);

	return &(botTree[newCenter]);

}


// Build a BST from the sorted extremes, return head of tree
bst_t* buildTree(const uint64_t *bot, const uint32_t m, const uint64_t tau, bst_t *botTree){

	bst_t *head; 
	uint64_t i;

	head = insertCenter(0, m-1, tau, bot, botTree);

	for (i=0; i<m; ++i){
		botTree[i].measNo = i;
	}

	return head;

}