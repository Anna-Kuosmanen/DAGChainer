/*
 *
 * SequenceGraph.cpp
 *
 * Created on: Sept 17th 2017
 *	Author: aekuosma
 *
 */

// TODO: Support for different lengths of pattern? Currently the pattern length has to be initialized when graph is read from files
// TODO: Exception checks, in case indexes > vector size are requested

#include "SequenceGraph.h"
#include "Tuple.h"

#include "Vertex.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>


void SequenceGraph::addToVertices(Vertex* v) {
	(this->vertices).push_back(v);
}

template<typename Out>
void SequenceGraph::split(const std::string &s, char delim, Out result) {
	std::stringstream ss;
	ss.str(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
	*(result++) = item;
	}
}

std::vector<std::string> SequenceGraph::split(const std::string &s, char delim) {
	std::vector<std::string> elems;
	split(s, delim, std::back_inserter(elems));
	return elems;
}


// Removes duplicates from the anchor list
// TODO Should this be done while searching for anchors or is it more effective here?
void SequenceGraph::removeDuplicates(std::vector<Tuple*> &anchors) {
	std::vector<Tuple*> tmp;

	for(int i=0;i<anchors.size();i++) {
		bool found = false;
		for(int j=0;j<tmp.size();j++) {
			if(anchors.at(i)->equals(*(tmp.at(j)))) {
				found = true;
				break;
			}
		}
 
		if(!found)
			tmp.push_back(anchors.at(i));
	}
	anchors = tmp;
}

SequenceGraph::SequenceGraph() {
}

// Need to define destructor to delete all those vertex
SequenceGraph::~SequenceGraph() {
	for(int i=0;i<this->getNoOfVertices();i++)
		delete vertices.at(i);
}

int SequenceGraph::getNoOfVertices() {
	return this->vertices.size();
}

// Returns i'th vertex
Vertex* SequenceGraph::getVertex(int i) {
	if(i< this->vertices.size())
		return this->vertices.at(i);
	else
		return NULL;
}

// Fills the i'th entry of array D (the length of the longest suffix of i'th prefix that matches the pattern)
// Also update the backpointers
// Returns true if label of v matched the pattern, false otherwise
// If i is outside pattern length, return false
bool SequenceGraph::fillDArray(Vertex* v, int i, std::string pattern) {
	int patlen = pattern.size();
	if(i >= patlen)
		return false;
	std::string prefix = pattern.substr(0,i+1);

	char ichar = prefix.at(prefix.size()-1);
	// the extension char doesn't match the label
	if(ichar != v->getLabel()) {
		return false;
	}
	//If checking the first prefix, it's of length 1, no check for extension
	else if(i==0) {
		v->updateDValue(i,1);
		return true;
	}
	std::vector<Vertex*> neighbors = v->getInNeighbors();

	int bestvalue = 0;
	Vertex* bestvertex = NULL;

	for(int j=0;j<neighbors.size();j++) {
		int compvalue = neighbors.at(j)->getDValue(i-1);

		// There should not be ties, otherwise the sequence graph would not have branched
		if(compvalue > bestvalue) {
			bestvalue = compvalue;
			bestvertex = neighbors.at(j);
		}
	}

	if(bestvertex != NULL) {
		v->addBackpointer(i, bestvertex);
	}

	v->updateDValue(i, bestvalue+1);

	return true;
		
}

// Reports matches that end at this vertex, for index i of D array
bool SequenceGraph::reportMatch(Vertex* v, int i, int threshold, std::vector<Tuple*> &results) {
	if(i >= 0 && v->getDValue(i) >= threshold) {
		std::vector<int> graphinterval;

		// Follow backpointers
		Vertex* current = v;
		int currenti = i;

		graphinterval.push_back(current->getId());

		// Where the path was cut
		int cuti = i;



		//TODO The following is clunky, make it nicer when have time
		// Check that the starting vertex has only one inneighbors
		if(current->getInNeighbors().size() > 1) {
			results.push_back(new Tuple(graphinterval, cuti,cuti));
			cuti = currenti-1;
			graphinterval.clear();
		}

		if(currenti > 0 && current->getBackpointer(currenti) != NULL) {
			do {	

				current = current->getBackpointer(currenti);
				currenti--;
				// If the graph branches, split at branchpoint (prevents overlaps)

				// More than one outneighbor, report existing path, excluding current vertex
				if(current->getOutNeighbors().size() > 1 && graphinterval.size() > 0) {
					std::vector<int> smallgraphinterval;
					// Reverse the path
					for(int k=graphinterval.size()-1;k>=0;k--)
						smallgraphinterval.push_back(graphinterval.at(k));
					results.push_back(new Tuple(smallgraphinterval, cuti-smallgraphinterval.size()+1,cuti));
					cuti = currenti;
					graphinterval.clear();

				}

				graphinterval.push_back(current->getId());

				// More than one inneighbors, report existing path, including this vertex
				if(current->getInNeighbors().size() > 1) {
					std::vector<int> smallgraphinterval;
					for(int k=graphinterval.size()-1;k>=0;k--) {
						smallgraphinterval.push_back(graphinterval.at(k));
					}
					results.push_back(new Tuple(smallgraphinterval, cuti-smallgraphinterval.size()+1,cuti));
					cuti = currenti-1;
					graphinterval.clear();
				}
			

		
			} while(currenti >= 0 && current->getBackpointer(currenti) != NULL);

		}

		// If there's anything that wasn't reported earlier
		if(graphinterval.size() > 0) {		

			// Reverse the graph interval vector, since it's from end to start now
			std::vector<int> rev;
			for(int k=graphinterval.size()-1;k>=0;k--)
				rev.push_back(graphinterval.at(k));

			results.push_back(new Tuple(rev, i-v->getDValue(i)+1,i-v->getDValue(i)+graphinterval.size()));
		}
		return true;
	}
	return false;
}

// Find all the anchors between a pattern and the sequence graph
void SequenceGraph::findAnchors(std::vector<Tuple*> &anchors, std::string pattern, int threshold) {

	// For every vertex
	for(int ver=0;ver<this->getNoOfVertices();ver++) {
		// For every length of prefix
		for(int i=0;i<=pattern.size();i++) {
			bool retval = this->fillDArray(this->vertices.at(ver),i, pattern);
			if(!retval) {
				// Label of v didn't match the pattern for this i
				// Report on every in-neighbor of v for index i-1, if D value is over 0
				std::vector<Vertex*> inneighbors = this->vertices.at(ver)->getInNeighbors();

				for(int j=0;j<inneighbors.size();j++) {
					if(inneighbors.at(j)->getDValue(i-1) > 0){
						bool retvalue = this->reportMatch(inneighbors.at(j),i-1, threshold, anchors);
					}
				}
			}
			// If v has no outneighbors, report on v (reportMatch will catch it if there's nothing to report)
			if(this->vertices.at(ver)->getOutNeighbors().size() == 0) {
				bool retvalue = this->reportMatch(this->vertices.at(ver) ,i, threshold, anchors);
			}		
		}
	}

	this->removeDuplicates(anchors);

}


// Graph file has first row the number of nodes (exons), then the edges (i+1'th row has the outgoing edges from i'th node)
// Nodes file has pairs of node number (first row) and the corresponding sequence (second row) for building the sequence graph (they need to be in order, labels are not checked currently!)
// Patlen is the length of the reads (needed to initialize the arrays in vertex objects)
void SequenceGraph::readSplicingGraphFile(std::string graphfile, std::string nodesfile, int patlen) {
	std::ifstream graphIn(graphfile.c_str());
	std::ifstream nodesIn(nodesfile.c_str());

	std::string line;
		
	int noOfExons = 0;

	// Read the first row to get the number of nodes (exons)
	if(graphIn.is_open()) {
		getline(graphIn,line);
		noOfExons = std::atoi(line.c_str());
	}

	std::vector<int> firstvertices;
	std::vector<int> lastvertices;

	int vertexCounter = 0;

	// Create the vertices. Save which vertex numbers are the first and last char of the splicing graph nodes
	if(nodesIn.is_open()) {
		for(int i=0;i<noOfExons;i++) {
			getline(nodesIn,line);
			std::string id = line;
			getline(nodesIn,line);
			std::string content = line;

			// One vertex per char
			// Add edge to previous if this isn't first vertex of some exon
			for(int j=0;j<content.size();j++) {
				Vertex* tmp = new Vertex(vertexCounter, content.at(j),patlen);
				this->addToVertices(tmp);
				if(j==0) {
					firstvertices.push_back(vertexCounter);
				}
				else {
					tmp->addEdgeFrom(this->vertices.at(this->vertices.size()-2));
				}
				if(j==content.size()-1) {
					lastvertices.push_back(vertexCounter);
				}
				vertexCounter++;
			}
		}
	}

	// Add the edges between exons
	// Delimiter is space
	if(graphIn.is_open()) {
		for(int i=0;i<noOfExons;i++) {
			getline(graphIn, line);
			if(line == "\n" || line == "")
				continue;
			std::vector<std::string> parts = split(line, ' ');
			std::vector<int> iparts;
			for(int j=0;j<parts.size();j++) {
				iparts.push_back(std::atoi(parts.at(j).c_str()));
			}

			for(int j=0;j<iparts.size();j++) {
				this->vertices.at(lastvertices.at(i))->addEdgeTo(this->vertices.at(firstvertices.at(iparts.at(j))));
			}
		}
	}

	graphIn.close();
	nodesIn.close();

}
