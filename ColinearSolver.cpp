/*
 * ColinearSolver.cpp
 *
 * Created on Oct 3rd 2017
 *	Author: Anna Kuosmanen
 *
 */

#include <iostream>
#include <fstream>

#include "ColinearSolver.h"
#include "RMaxQTree.h"

#include <lemon/lgf_reader.h>
// These are included through the MPC header files
//#include <lemon/list_graph.h>
//#include <lemon/adaptors.h>
//#include <lemon/connectivity.h>
//#include <lemon/lgf_writer.h>
//#include <climits> // INT_MAX

#include "MC-MPC/util/utils.h"
#include "MC-MPC/decomposer/MPC.h"
#include "MC-MPC/decomposer/decomposition.h"

// For sorting the key pairs
struct less_than_key {
	inline bool operator() (const std::pair<int,int> key1, std::pair<int,int> key2) {
		if(key1.first != key2.first)
			return (key1.first < key2.first);
		else
			return (key1.second < key2.second);

	}
};

int ColinearSolver::getPathCoverSize() {
	return this->pathcover.size();
}

// Converts the SequenceGraph objects into flow graph and writes it to filename
void ColinearSolver::convertSequenceGraphToFlowGraph(std::string filename) {

	std::ofstream out(filename.c_str());

	out << "@nodes" << std::endl << "label" << std::endl;

	for(int i=0;i<this->SGraph->getNoOfVertices();i++) {
		out << i << std::endl;
	}

	out << "@arcs" << std::endl;

	out << "\t\tlabel\tweight" << std::endl;

	int arccount = 1;
	for(int i=0;i<this->SGraph->getNoOfVertices();i++) {

		std::vector<Vertex*> outneighbors = this->SGraph->getVertex(i)->getOutNeighbors();

		for(unsigned j=0;j<outneighbors.size();j++) {

			// From this vertex (id is the number of the vertex)...
			out << i << " ";

			//...to these vertices
			out << outneighbors.at(j)->getId() << " ";

			// And arc ID
			out << arccount << " ";
			arccount++;

			// and weight (1)
			out << "1" << std::endl;

		}
	}

}

// Read the flow graph from filename and solve for paths, then convert to paths in the original graph
// (Drop source and sink and -1 on all indexes)
void ColinearSolver::solveForPaths(std::string filename) {

	// start of Topi's code
	ListDigraph graph;

	ListDigraph::NodeMap<int> node_labels(graph);
	ListDigraph::ArcMap<int> arc_labels(graph);
	ListDigraph::ArcMap<int> arc_weights(graph);

	digraphReader(graph, filename)
		.nodeMap("label", node_labels)
		.arcMap("label", arc_labels)
		.arcMap("weight", arc_weights)
		.run();


	ListDigraph::Node s = add_source(graph);
	ListDigraph::Node t = add_sink(graph);
	
	ListDigraph::ArcMap<int> minflow(graph);
	find_minflow(graph, minflow, s, t);

	int num_paths = 0;
	for(ListDigraph::OutArcIt o(graph, s); o != INVALID; ++o){
		num_paths += minflow[o];
	}
	//Extract paths from the graph according to the minFlow

	vector<ListDigraph::Node>* paths = new vector<ListDigraph::Node>[num_paths];
//	ListDigraph::NodeMap<int*> reachable(graph);

	//paths are picked up one path at time
	for (int i = 0; i < num_paths; ++i){
		ListDigraph::Node node = s;
		int path_index = 0;
		while(node != t){
			for (ListDigraph::OutArcIt o(graph, node); o != INVALID; ++o){
				if(minflow[o] > 0){
					minflow[o] -= 1;
					paths[i].push_back(node);
					node = graph.target(o);
					break;
				}
			}
			path_index++;
		}

		paths[i].push_back(t);
	}
	// End of Topi's code

	// Convert to paths and save to "paths" variable
	for(int k=0;k<num_paths;k++) {
		std::vector<int> path;
		// Drop source and sink nodes
		for(unsigned l=1;l<paths[k].size()-1;l++) {
			// Map back to sequencegraph vertices
			path.push_back(this->SGraph->getVertex(graph.id(paths[k].at(l)))->getId());
		}
		this->pathcover.push_back(path);
	}

	// Clean-up
	delete [] paths;
}


// Computes forward links for all vertices
void ColinearSolver::computeForward() {

	int index[this->SGraph->getNoOfVertices()][this->pathcover.size()];
	int distlast2reach[this->SGraph->getNoOfVertices()][this->pathcover.size()];

	// Init the paths and forward vector
	for(int i=0;i<this->SGraph->getNoOfVertices();i++) {
		std::vector<int> tmp;
		this->pathsforv.push_back(tmp);
		std::vector<std::pair<int,int> > tmp2;
		this->forward.push_back(tmp2);
	}

	for(unsigned i=0;i<this->pathcover.size();i++) {
		vector<int> path = pathcover.at(i);
		for(unsigned j=0;j<path.size();j++) {
			int v = path.at(j);
			this->pathsforv[v].push_back(i);
			index[v][i] = j;
		}
	}


	int last2reach[this->SGraph->getNoOfVertices()][this->pathcover.size()];

	// Init all to -1
	for(int v=0;v<this->SGraph->getNoOfVertices();v++) {
		for(unsigned i=0;i<pathcover.size();i++) {
			last2reach[v][i] = -1;
			distlast2reach[v][i] = -1;
		}
	}

	for(int v=0;v<this->SGraph->getNoOfVertices();v++) {
		std::vector<int> vpaths = pathsforv[v];
		for(unsigned i=0;i<pathcover.size();i++) {
			// i in paths[v]
			if(find(vpaths.begin(),vpaths.end(),i)!=vpaths.end()) {
				last2reach[v][i] = index[v][i];
				distlast2reach[v][i] = 0;
			}
			// i not in paths[v]
			else {
				vector<Vertex*> inneighbors = SGraph->getVertex(v)->getInNeighbors();
				for(unsigned k=0;k<inneighbors.size();k++) {
					int u = inneighbors.at(k)->getId();
					if(last2reach[u][i] > last2reach[v][i]) {
						last2reach[v][i] = last2reach[u][i];
						distlast2reach[v][i] = distlast2reach[u][i]+1;
					}
				}
			
			}
		}
	}

	for(int v=0;v<this->SGraph->getNoOfVertices();v++) {
		for(unsigned i=0;i<pathcover.size();i++) {
			if(last2reach[v][i] != -1) {
				if(distlast2reach[v][i] <= this->maxdistance) {
					forward.at(pathcover.at(i).at(last2reach[v][i])).push_back(std::make_pair(v,i));
				}
			}
		}
	}
}

// +1 every key coordinate to shift from 0-based to 1-based (otherwise the key 0 causes issues)
ColinearChain ColinearSolver::solveForAnchors(std::vector<Tuple*> &M) {

	// If there are no anchors, return failure
	if(M.size() == 0)
		return ColinearChain();

	RMaxQTree* ITrees = new RMaxQTree[this->pathcover.size()];
	RMaxQTree* TTrees = new RMaxQTree[this->pathcover.size()];
		
	// Get the key list from anchors: all pairs (M[j].d,j)
	std::vector<std::pair<int,int> > keys;
	keys.push_back(std::make_pair(0,-1));
		
	for(unsigned i=0;i<M.size();i++) {
		keys.push_back(std::make_pair(M.at(i)->d+1,i));
	}
	
	std::sort(keys.begin(),keys.end(), less_than_key());
	
	// Change into array for the RMaxQTree
	std::pair<int,int> keysarray[keys.size()];
	
	for(unsigned i=0;i<keys.size();i++)
		keysarray[i] = keys.at(i);
		
	// Init the trees
	for(unsigned i=0;i<this->pathcover.size();i++) {
		ITrees[i].fillRMaxQTree(keysarray, keys.size());
		TTrees[i].fillRMaxQTree(keysarray, keys.size());
		ITrees[i].update(std::make_pair(0,-1),-1,0);
		TTrees[i].update(std::make_pair(0,-1),-1,0);
	}

	std::vector<int>* start = new std::vector<int>[this->SGraph->getNoOfVertices()];
	std::vector<int>* end = new std::vector<int>[this->SGraph->getNoOfVertices()];
		
	for(unsigned j=0;j<M.size();j++) {
		start[M.at(j)->PFirst].push_back(j);
		end[M.at(j)->PLast].push_back(j);
	}
	
	
	for(int v=0;v<this->SGraph->getNoOfVertices();v++) {

		// Remove the coverage of tuples that are too far (i.e. set their best coverage to -infinity)
		std::vector<int> pathsv = this->pathsforv[v];

		for(unsigned i=0;i<pathsv.size();i++) {
			// Find how many'th v is on this path
			int indexOnPath = -1;
			for(unsigned j=0;j<pathcover.at(pathsv.at(i)).size();j++) {
				if(pathcover.at(pathsv.at(i)).at(j) == v) {
					indexOnPath = j;
					break;
				}
			}
			// Nothing is too far back
			if(indexOnPath < this->maxdistance)
				continue;

			std::vector<int> deltuples = end[pathcover.at(pathsv.at(i)).at(indexOnPath-this->maxdistance)];

			for(unsigned t=0;t<deltuples.size();t++) {
				TTrees[pathsv.at(i)].update(std::make_pair(M.at(deltuples.at(t))->d+1,deltuples.at(t)),-1,negative_infinity);
				ITrees[pathsv.at(i)].update(std::make_pair(M.at(deltuples.at(t))->d+1,deltuples.at(t)),-1,negative_infinity);		
			}

		}


		std::vector<std::pair<int,int> > forwardv = this->forward[v];

		// First process possible forward links to self
		// This is strictly to handle cases where one path is prefix of another and has length 1 (as the loop below has to update coverage or path of length 1 is ignored)
		for(unsigned f=0;f<forwardv.size();f++) {
			std::vector<int> startw = start[forwardv.at(f).first];
			
			for(unsigned j=0;j<startw.size();j++) {

				if(M.at(startw.at(j))->PFirst != v)
					continue;

				std::pair<int,int> atuple = TTrees[forwardv.at(f).second].query(0,M.at(startw.at(j))->c+1-1);
				std::pair<int,int> btuple = ITrees[forwardv.at(f).second].query(M.at(startw.at(j))->c+1, M.at(startw.at(j))->d+1);

				int Ca = (M[startw.at(j)]->d-M[startw.at(j)]->c+1) + atuple.second;
				int Cb = M[startw.at(j)]->d + btuple.second;
					
				if(Ca > M.at(startw.at(j))->C && Ca >= Cb) {
					M.at(startw.at(j))->C = Ca;
					if(atuple.first != -1)
						M.at(startw.at(j))->previous = M.at(atuple.first);
				}
				else if(Cb > M.at(startw.at(j))->C) {
					M.at(startw.at(j))->C = Cb;
					if(btuple.first != -1)
						M.at(startw.at(j))->previous = M.at(btuple.first);
				}
			}
		}


		// For every tuple whose path ends at this vertex, update the trees for all the paths on which this vertex lies
		for(unsigned j=0;j<end[v].size();j++) {

//			std::vector<int> pathsv = this->pathsforv[v];
			for(unsigned i=0;i<pathsv.size();i++) {

				TTrees[pathsv.at(i)].update(std::make_pair(M.at(end[v].at(j))->d+1,end[v].at(j)), end[v].at(j),M.at(end[v].at(j))->C);
				ITrees[pathsv.at(i)].update(std::make_pair(M.at(end[v].at(j))->d+1,end[v].at(j)), end[v].at(j),M.at(end[v].at(j))->C-(M.at(end[v].at(j))->d));
			}
		}
		
		// And then follow the forward links to update the coverages for tuples whose first node is in the forward set
		for(unsigned f=0;f<forwardv.size();f++) {
			std::vector<int> startw = start[forwardv.at(f).first];
			
			for(unsigned j=0;j<startw.size();j++) {

				// Skip if it's link to self, already processed
				if(M.at(startw.at(j))->PFirst == v)
					continue;

				std::pair<int,int> atuple = TTrees[forwardv.at(f).second].query(0,M.at(startw.at(j))->c+1-1);
				std::pair<int,int> btuple = ITrees[forwardv.at(f).second].query(M.at(startw.at(j))->c+1, M.at(startw.at(j))->d+1);

				int Ca = (M[startw.at(j)]->d-M[startw.at(j)]->c+1) + atuple.second;
				int Cb = M[startw.at(j)]->d + btuple.second;
					
				if(Ca > M.at(startw.at(j))->C && Ca >= Cb) {
					M.at(startw.at(j))->C = Ca;
					if(atuple.first != -1)
						M.at(startw.at(j))->previous = M.at(atuple.first);
				}
				else if(Cb > M.at(startw.at(j))->C) {
					M.at(startw.at(j))->C = Cb;
					if(btuple.first != -1)
						M.at(startw.at(j))->previous = M.at(btuple.first);
				}
			}
		}
	}
	
	std::vector<Tuple*> tempsolution;
	
	// Find max C[j]
	int maxvalue = 0;
	int maxindex = -1;
	for(unsigned i=0;i<M.size();i++) {
		if(M.at(i)->C > maxvalue) {
			maxvalue = M.at(i)->C;
			maxindex = i;
		}
	}
	// ... and backtrack
	Tuple* currentTuple = M.at(maxindex);
	tempsolution.push_back(currentTuple);
	
	while(currentTuple->previous != NULL) {
		currentTuple = currentTuple->previous;
		tempsolution.push_back(currentTuple);
	}
	
	std::vector<Tuple*> solution;

	// It's backwards, so reverse
	for(int i=tempsolution.size()-1;i>=0;i--)
		solution.push_back(tempsolution.at(i));
	

	// Clean-up
	delete [] start;
	delete [] end;
	delete [] TTrees;
	delete [] ITrees;

	if(solution.size() > 0)
		return ColinearChain(solution, maxvalue);
	else
		return ColinearChain();
	
}

// TODO Add check that SequenceGraph isn't empty (first giving the graph to solver and then reading it doesn't work)
// TODO Relative path for the flowgraph.tmp, it's not compatible with running multiple processes in one directory
ColinearSolver::ColinearSolver(SequenceGraph* &SGraph, int maxdistance) {
	this->SGraph = SGraph;
	this->maxdistance = maxdistance;

	// Compute the path cover and forward links for this SequenceGraph
	std::string tempgraphfile = "flowgraph.tmp";

	this->convertSequenceGraphToFlowGraph(tempgraphfile);
	this->solveForPaths(tempgraphfile);

	this->computeForward();

}

ColinearSolver::~ColinearSolver() {
	this->clearAnchors();
	this->lastanchors.clear();
	pathcover.clear();
	pathsforv.clear();
	forward.clear();
}

void ColinearSolver::clearAnchors() {
	for(unsigned i=0;i<this->lastanchors.size();i++)
		delete this->lastanchors.at(i);
}


// Solves the co-linear chaining problem for the given read with seed threshold of "minthres" ( = minimum MEM length)
// If the last parameter is true, check both read and its reverse complement, otherwise just the read
ColinearChain ColinearSolver::solve(FastaEntry read, std::string gcsa_file, std::string lcp_file, int minthres, bool bothways) {
	this->clearAnchors();

	std::string id = read.id;
	std::vector<Tuple*> anchors;
	std::vector<Tuple*> revanchors;

	// Check both sequence and its reverse complement, and pick the one with better score
	if(bothways) {

		// TODO Would be nice to create the index objects here and pass them on, but for some reason that causes segmentation fault in gcsa->alpha.char2comp (the object exists, but is not initialized?)
		SGraph->findAnchorsGCSA2(anchors, revanchors, gcsa_file, lcp_file, read.seq, minthres); 

		ColinearChain chain = this->solveForAnchors(anchors);
		ColinearChain revchain = this->solveForAnchors(revanchors);

		if(chain.coverageScore >= revchain.coverageScore) {

			lastanchors = anchors;
			for(unsigned i=0;i<revanchors.size();i++)
				delete revanchors.at(i);
		
			return chain;
		}
		else {
			lastanchors = revanchors;
			for(unsigned i=0;i<anchors.size();i++)
				delete anchors.at(i);

			return revchain;
		}
	}
	else {
		SGraph->findAnchorsGCSA2(anchors, revanchors, gcsa_file, lcp_file, read.seq, minthres); 

		ColinearChain chain = this->solveForAnchors(anchors);

		for(unsigned i=0;i<revanchors.size();i++)
			delete revanchors.at(i);

		lastanchors = anchors;

		return chain;
	}
}
