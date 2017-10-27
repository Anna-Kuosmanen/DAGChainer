/*
 *
 * ColinearSolver.h
 *
 * Created on Oct 3rd 2017
 *	Author: aekuosma
 *
 */

#ifndef COLINEARSOLVER_H_
#define COLINEARSOLVER_H_


#include "SequenceGraph.h"
#include "Vertex.h"

#include <vector>

class ColinearSolver {

private:

	SequenceGraph* SGraph;
	std::vector<std::string> patterns;

	// Contains the vertex sequence for each path
	std::vector<std::vector<int> > pathcover;

	// Contains, for each vertex v, the paths in lies on
	std::vector<std::vector<int> > pathsforv;

	// forward[u] from the paper
	std::vector<std::vector<std::pair<int, int> > > forward;


public:

	int getPathCoverSize();

	// Converts the SequenceGraph objects into flow graph and writes it to filename
	void convertSequenceGraphToFlowGraph(std::string filename);
	
	// Read the flow graph from filename and solve for paths, then convert to paths in the original graph
	void solveForPaths(std::string filename);

	// Computes forward links for all vertices
	void computeForward();
	
	std::vector<Tuple*> solveForAnchors(std::vector<Tuple*> M);

	ColinearSolver(SequenceGraph* SGraph);
	
	// Solves the co-linear chaining problem
	void solve(std::string patternfile, std::string outputfile, int threshold);

};

#endif /* COLINEARSOLVER_H_ */
