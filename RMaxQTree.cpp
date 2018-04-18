/*
 * RMaxQTree.cpp
 *
 * Created Oct 7th 2017
 *	Author: Anna Kuosmanen
 *
 */

#include "RMaxQTree.h"

#include <vector>
#include <limits>
#include <math.h>
#include <iostream>

//int negative_infinity = - std::numeric_limits<int>::infinity();
int negative_infinity = -10000000;

void RMaxQTree::init(int node, int b, int e, int *keys) {
	// leaf
	if (b == e)
		tree[node].setValues(keys[b],-1, negative_infinity);
	else {
		// split
		init(2 * node, b, (b + e) / 2, keys);
		init(2 * node + 1, (b + e) / 2 + 1, e, keys);
		// propagate up
		if (tree[2 * node].Cj >= tree[2 * node + 1].Cj)
			tree[node] = tree[2 * node];
		else
			tree[node] = tree[2 * node + 1];
	}
}
	
	
// Called by public function update with only two parameters
void RMaxQTree::updateTree(int key, int j, int Cj, int node, int b, int e) {
	// Ended in a leaf
	if (b == e) {
		if(tree[node].key == key) {
			tree[node].j = j;
			tree[node].Cj = Cj;
		}
	}
	// Search
	else {	
		int mid = (b + e) / 2;
		if (key <= keys[mid])
			updateTree(key, j, Cj, 2 * node, b, mid);
		else
			updateTree(key, j, Cj, 2 * node + 1, mid + 1, e);
		// And propagate back up
		if (tree[2 * node].Cj >= tree[2 * node + 1].Cj)
			tree[node] = tree[2 * node];
		else
			tree[node] = tree[2 * node + 1];
	}
}

// Called by public function query with only two parameters
// Returns pair (j,C[j])
// Note here that i and j are values of _keys_, whereas b and e are indexes in the keys array!
std::pair<int,int> RMaxQTree::queryTree(int i, int j, int node, int b, int e) {
	// bad interval
	if (i > keys[e] || j < keys[b])
		return std::make_pair(-1,negative_infinity);
	// good interval
	if (keys[b] >= i && keys[e] <= j)
		return std::make_pair(tree[node].j,tree[node].Cj);
	// check left and right subtree
	std::pair<int,int> left = queryTree(i, j, 2 * node, b, (b + e) / 2);
	std::pair<int,int> right = queryTree(i, j, 2 * node + 1, (b + e) / 2 + 1, e);
	if (left.second == negative_infinity)
		return right;
	if (right.second == negative_infinity)
		return left;
	if (left.second >= right.second)
		return left;
	return right;
}
	

// Empty constructor for creating arrays
RMaxQTree::RMaxQTree() {}

RMaxQTree::~RMaxQTree() {
	delete [] this->tree;
}

// For filling the empty RMaxQTrees
void RMaxQTree::fillRMaxQTree(int *keys, int keyLen) {
	this->keyLen = keyLen;
	this->keys = keys;
	this->treeLen = 2 << (int)ceil(log2(keyLen));
	this->tree = new TreeNode[treeLen];
	init(1, 0, keyLen - 1, keys);	
}	

RMaxQTree::RMaxQTree(int *keys, int keyLen) {
	this->keyLen = keyLen;
	this->keys = keys;
	this->treeLen = 2 << (int)ceil(log2(keyLen));
	this->tree = new TreeNode[treeLen];
	init(1, 0, keyLen - 1, keys);
}

void RMaxQTree::update(int key, int j, int Cj) {
	// Start with node 1, go over the whole tree till find the key
	this->updateTree(key, j, Cj, 1, 0, keyLen - 1);
}
	
std::pair<int,int> RMaxQTree::query(int start, int end) {
	// Go over the whole tree till find the query
	return this->queryTree(start, end, 1, 0, keyLen - 1);
}
 
