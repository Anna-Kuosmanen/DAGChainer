/*
 * ColinearSolver.h
 *
 * Created on Oct 3rd 2017
 *	Author: Anna Kuosmanen
 *
 */

#ifndef COLINEARSOLVER_H_
#define COLINEARSOLVER_H_

#include "SequenceGraph.h"
#include "Vertex.h"
#include "ColinearChain.h"
#include "FastaEntry.h"
#include "MaximalExactMatch.h"

#include <vector>

class ColinearSolver {

private:

	SequenceGraph* SGraph;

	// Contains the vertex sequence for each path
	std::vector<std::vector<int> > pathcover;

	// Contains, for each vertex v, the paths it lies on
	std::vector<std::vector<int> > pathsforv;

	// forward[u] from the RECOMB2018 paper
	std::vector<std::vector<std::pair<int, int> > > forward;

	// For clean-up purposes, otherwise the old anchor pointers are lost
	std::vector<Tuple*> lastanchors;

	// For using GCSA2 index
	gcsa::GCSA gcsa;
	gcsa::LCPArray lcp;

public:

	int getPathCoverSize();

	// Converts the SequenceGraph objects into flow graph and writes it to filename
	void convertSequenceGraphToFlowGraph(std::string filename);
	
	// Read the flow graph from filename and solve for paths, then convert to paths in the original graph
	void solveForPaths(std::string filename);

	// Computes forward links for all vertices
	void computeForward();
	
	// Finds the best chain for set of anchors M
	ColinearChain solveForAnchors(std::vector<Tuple*> &M);

	// Read the SequenceGraph from splicing graph file first,
	// reading after solver is created doesn't work
	ColinearSolver(SequenceGraph* &SGraph);
	
	~ColinearSolver();

	// Cleans up
	void clearAnchors();

	// Solves the co-linear chaining problem for the given read
	// If last parameter is true, checks both the read and its reverse complement
	// Otherwise check only the read
	ColinearChain solve(FastaEntry read, std::string gcsa_file, std::string lcp_file, int minthres, bool bothways);

};

#endif /* COLINEARSOLVER_H_ */
