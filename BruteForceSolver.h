/*
 * BruteForceSolver.h
 *
 * Created on Oct 23rd 2017
 *	Author: aekuosma
 *
 */

#ifndef _BRUTEFORCESOLVER_H_
#define _BRUTEFORCESOLVER_H_

#include "SequenceGraph.h"

class BruteForceSolver {

private:

	SequenceGraph* SGraph;


public:

	BruteForceSolver(SequenceGraph* &SGraph);

	// The recursive function that does DFS and updates Tuple orig's coverage as it goes
	void DFSUpdate(Tuple* orig, Vertex* v, std::vector<bool> &covered, std::vector<std::vector<int> > &end, std::vector<Tuple*> &M);	

	void solveForAnchors(std::vector<Tuple*> &M, std::vector<Tuple*> &solution);

	// Threshold is the minimum length of the anchors
	void solve(std::string patternfile, std::string outputfile, int threshold);

};
#endif // _BRUTEFORCESOLVER_H_
