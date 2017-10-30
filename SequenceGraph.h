/*
 *
 * SequenceGraph.h
 *
 * Created on: Sept 17th 2017
 *	Author: aekuosma
 *
 */

#ifndef SEQUENCEGRAPH_H_
#define SEQUENCEGRAPH_H_

#include <vector>
#include <iostream>

#include "Vertex.h"
#include "Tuple.h"


class SequenceGraph {

public:
	std::vector<Vertex*> vertices; // Vertices in topological order. 
				  // Currently the order comes from the file format (splicing graph files are ordered by genomic coordinates).

	// Actual creation is done in the reading function
	SequenceGraph();

	void addToVertices(Vertex* v);

	template<typename Out>
	void split(const std::string &s, char delim, Out result);

	std::vector<std::string> split(const std::string &s, char delim);


	// Removes duplicates from the anchor list
	// TODO Should this be done while searching for anchors or is it more effective here?
	void removeDuplicates(std::vector<Tuple*> &anchors);

	// Destructor
	~SequenceGraph();
	
	int getNoOfVertices();

	// Returns i'th vertex
	Vertex* getVertex(int i);

	// Fills the i'th entry of array D (the length of the longest suffix of i'th prefix that matches the pattern)
	// Also update the backpointers
	// Returns true if label of v matched the pattern, false otherwise
	// If i is outside pattern length, return false
	bool fillDArray(Vertex* v, int i, std::string pattern);

	void clearDArray(int patlen);

	// Reports matches that end at this vertex, for index i of D array
	// Report is in the form of ((start in pattern, end in pattern),(list of path nodes in graph))
		
	bool reportMatch(Vertex* v, int i, int threshold, std::vector<Tuple*> &results);

	// Find and return all the anchors between a pattern and the sequence graph
		
	// Note that there's no reporting if the length of the prefix is 1, it has to be handled here
	void findAnchors(std::vector<Tuple*> &results, std::string pattern, int threshold);

	// Graph file has first row the number of nodes (exons), then the edges (i+1'th row has the outgoing edges from i'th node)
	// Nodes file has pairs of node number (first row) and the corresponding sequence (second row) for building the sequence graph (they need to be in order, labels are not checked currently!)
	// Patlen is the length of the reads (needed to initialize the arrays in vertex objects)
	void readSplicingGraphFile(std::string graphfile, std::string nodesfile, int patlen);
	

};

#endif /* SEQUENCEGRAPH_H_ */
