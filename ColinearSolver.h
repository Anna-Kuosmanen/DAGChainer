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

	// Converts the SequenceGraph objects into flow graph and writes it to filename
	void convertSequenceGraphToFlowGraph(std::string filename);
	
	// Read the flow graph from filename and solve for paths, then convert to paths in the original graph
	void solveForPaths(std::string filename);

	// Computes forward links for all vertices
	void computeForward();
	
	std::vector<Tuple*> colinearChain(std::vector<Tuple*> anchors);


public:

	// Note: Patterns have to be the same length
	ColinearSolver(std::string graphFile, std::string nodesFile, int patlen);
	
	// Solves the co-linear chaining problem (patterns are read from the file, with one pattern per row)
	// Outputs to the given file with first row having the pattern, following row the number of anchors,
	// and then all the anchor tuples, one per row ((patstart,patend),(path in graph))
	void solve(std::string patternfile, std::string outputfile, int threshold);

};

#endif /* COLINEARSOLVER_H_ */
