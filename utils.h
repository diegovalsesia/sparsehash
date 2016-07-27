#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef struct bst{

	uint64_t measNo;
	uint64_t botVal;
	uint64_t topVal;
	struct bst *leftPtr;
	struct bst *rightPtr;

} bst_t;


// Build a BST from the sorted extremes, return head of tree
bst_t* buildTree(const uint64_t *bot, const uint32_t m, const uint64_t tau, bst_t *botTree);