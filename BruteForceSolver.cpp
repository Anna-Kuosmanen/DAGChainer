/*
 * BruteForceSolver.cpp
 *
 * Created on Oct 23rd 2017
 *	Author: aekuosma
 *
 */

#include "BruteForceSolver.h"

#include <iostream>
#include <fstream>
#include <algorithm>


BruteForceSolver::BruteForceSolver(SequenceGraph* &SGraph) {
	this->SGraph = SGraph;
}


// Orig is the original tuple for which search for the best coverage
// v is the current vertex under processing
// covered lists which nodes have been covered
// end lists all the tuples that end at vertex at position i
// M is the list of anchors
void BruteForceSolver::DFSUpdate(Tuple* orig, Vertex* v, std::vector<bool> &covered, std::vector<std::vector<int> > &end, std::vector<Tuple*> &M) {
	covered.at(v->getId()) = true;
	// Check if some tuple ends at this vertex
	for(int i=0;i<end.at(v->getId()).size();i++) {
		Tuple* candtuple = M.at(end.at(v->getId()).at(i));
		int C = -1;
		if(candtuple->d < orig->c)
			C = candtuple->C + (orig->d-orig->c+1);
		else if(orig->c <=candtuple->d && candtuple->d <= orig->d)
			C = candtuple->C +(orig->d-candtuple->d);

		// If this had better coverage, update originals and the pointer
		if(C > orig->C) {
			orig->C = C;
			orig->previous = candtuple;
		}
	}

	// Recur for every outneighbor in the reverse graph, i.e. for every inneighbor
	std::vector<Vertex*> outneighbors = v->getInNeighbors();

	for(int i=0;i<outneighbors.size();i++) {
		int othercov = 0;
		if(!(covered.at(outneighbors.at(i)->getId())))
			DFSUpdate(orig, outneighbors.at(i), covered, end, M);

	}

}

// Applies the brute-force algorithm on the given anchors
void BruteForceSolver::solveForAnchors(std::vector<Tuple*> &M, std::vector<Tuple*> &solution) {

	std::vector<std::vector<int> > end;

	// For depth-first search
	std::vector<bool> covered;

	for(int i=0;i<this->SGraph->getNoOfVertices();i++) {
		std::vector<int> tmp;
		end.push_back(tmp);
		covered.push_back(false);
	}

	// At end[i] add the indexes of all pairs whose path ends at i
	for(int j=0;j<M.size();j++) {
		end.at(M.at(j)->PLast).push_back(j);
	}

	// Sort
	std::vector<Tuple*> temp;

	for(int i=0;i<this->SGraph->getNoOfVertices();i++) {
		std::vector<int> ends_for_one = end.at(i);
		for(int j=0;j<ends_for_one.size();j++) {
			temp.push_back(M.at(ends_for_one.at(j)));
		}

	}

	M = temp;

	// Re-populate end
	for(int i=0;i<this->SGraph->getNoOfVertices();i++)
		end.at(i).clear();

	for(int j=0;j<M.size();j++) {
		end.at(M.at(j)->PLast).push_back(j);
	}


	for(int j=0;j<M.size();j++) {

		// The coverage if this is the first of the chain
		M.at(j)->C = M.at(j)->d - M.at(j)->c +1;

		Vertex* w = this->SGraph->getVertex(M.at(j)->PFirst);
		std::vector<Vertex*> neighbors = w->getInNeighbors();

		// Depth-first search for every inneighbor (outneighbor in reverse) of startpoint of this tuple, updates directly to M.at(j)->C
		for(int k=0;k<neighbors.size();k++) {

			DFSUpdate(M.at(j), neighbors.at(k), covered, end, M);
		}
		//Set covered to false again for nex tuple
		for(int i=0;i<this->SGraph->getNoOfVertices();i++)
			covered.at(i) = false;
	}

	std::vector<Tuple*> tempsolution;

	// Find max C[j]
	int maxvalue = 0;
	int maxindex = -1;
	for(int i=0;i<M.size();i++) {
		if(M.at(i)->C > maxvalue) {
			maxvalue = M.at(i)->C;
			maxindex = i;
		}
	}
	// ... and backtrack
	Tuple* currentTuple = M.at(maxindex);
	tempsolution.push_back(currentTuple);

	while(currentTuple->previous != NULL) {
		currentTuple =currentTuple->previous;
		tempsolution.push_back(currentTuple);
	}

	// It's backwards, so reverse

	for(int i=solution.size()-1;i>=0;i--)
		solution.push_back(tempsolution.at(i));

}



void BruteForceSolver::solve(std::string patternfile, std::string outputfile, int threshold){

        std::ifstream patin(patternfile.c_str());
        std::ofstream out(outputfile.c_str());

	std::string line;

        while(getline(patin, line)) {
                std::vector<Tuple*> anchors;
		SGraph->findAnchors(anchors, line, threshold);
                std::vector<Tuple*> chain;
		solveForAnchors(anchors, chain);

                for(int i=0;i<chain.size();i++)
                        out << chain.at(i)->toString() << ", ";

                out << std::endl;

		// Clean-up
		for(int i=0;i<anchors.size();i++)
			delete anchors.at(i);
        }
}
