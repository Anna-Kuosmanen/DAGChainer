/*
 * RMaxQTree.h
 *
 * Created Oct 7th 2017
 *	Author: Anna Kuosmanen
 *
 * A Range Maximum Query Tree.
 *
 * A static implementation. All the keys must be given at the construction.
 */

#ifndef RMaxQTREE_H_
#define RMaxQTREE_H_

#include <utility>

//static int negative_infinity = - std::numeric_limits<int>::infinity();
static int negative_infinity = -10000000;


class TreeNode {
public:
	std::pair<int,int> key; // keys are pairs (M[j].d,j)
	int Cj; // The best coverage
	int j; // The Tuple index that the best coverage corresponds to

	// Dummy needed for initializing the array
	TreeNode() {
	}
	
	TreeNode(std::pair<int,int> key, int j, int Cj) {
		this->key = key;
		this->Cj = Cj;
		this->j = j;
	}

	void setValues(std::pair<int,int> key, int j, int Cj) {
		this->key = key;
		this->Cj = Cj;
		this->j = j;
	}
};

class RMaxQTree {

private:
	TreeNode *tree;
	std::pair<int,int> *keys;
	int treeLen, keyLen;

	void init(int node, int b, int e, std::pair<int,int> *keys);

	// Recursive helper function
	void updateTree(std::pair<int,int> key, int j, int Cj, int node, int b, int e);

	// Recursive helper function
	std::pair<int,int> queryTree(int i, int j, int node, int b, int e);

public:

	// Empty constructor for creating arrays
	RMaxQTree();

	~RMaxQTree();

	// For filling the empty RMaxQTrees
	void fillRMaxQTree(std::pair<int,int> *keys, int keyLen);

	RMaxQTree(std::pair<int,int> *keys, int keyLen);

	// Update the node with key "key" with the given values
	void update(std::pair<int,int> key, int j, int Cj);

	// Return pair (j,C[j])
	std::pair<int,int> query(int start, int end);
};

#endif /* _RMaxQTREE_H_ */
