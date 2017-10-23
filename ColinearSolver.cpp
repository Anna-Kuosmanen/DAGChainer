/*
 *
 * ColinearSolver.cpp
 *
 * Created on Oct 3rd 2017
 *	Author: aekuosma
 *
 */

#include "ColinearSolver.h"
 
#include "SequenceGraph.h"
#include "RMaxQTree.h"
//#include "Tuple.h"

#include <lemon/list_graph.h>
#include <lemon/adaptors.h>
#include <lemon/connectivity.h>
#include <lemon/lgf_reader.h>
#include <lemon/lgf_writer.h>
#include <climits> // INT_MAX

#include "MC-MPC/util/utils.h"
#include "MC-MPC/decomposer/MPC.h"
#include "MC-MPC/decomposer/decomposition.h"


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

		for(int j=0;j<outneighbors.size();j++) {

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
	ListDigraph::NodeMap<int*> reachable(graph);

	//paths are picked up one path at time
	for (int i = 0; i < num_paths; ++i){
		ListDigraph::Node node = s;
		int path_index = 0;
		while(node != t){
			if(reachable[node] == NULL){
				reachable[node] = (int*) calloc(num_paths, sizeof(int));
				for (int i = 0; i < num_paths; ++i){
					reachable[node][i] = INT_MAX;
				}
			}
			reachable[node][i] = path_index;

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

		if(reachable[t] == NULL){
			reachable[t] = (int*) calloc(num_paths, sizeof(int));
		}
		reachable[t][i] = path_index;
		paths[i].push_back(t);
	}
	// End of Topi's code

	// Convert to paths and save to "paths" variable
	for(int k=0;k<num_paths;k++) {
		std::vector<int> path;
		// Drop source and sink nodes
		for(int l=1;l<paths[k].size()-1;l++) {
			// Map back to sequencegraph vertices
			path.push_back(this->SGraph->getVertex(graph.id(paths[k].at(l)))->getId());
		}
		this->pathcover.push_back(path);
	}
}


// Computes forward links for all vertices
void ColinearSolver::computeForward() {

	int index[this->SGraph->getNoOfVertices()][this->pathcover.size()];

	// Init the paths and forward vector
	for(int i=0;i<this->SGraph->getNoOfVertices();i++) {
		std::vector<int> tmp;
		this->pathsforv.push_back(tmp);
		std::vector<std::pair<int,int> > tmp2;
		this->forward.push_back(tmp2);
	}

	for(int i=0;i<this->pathcover.size();i++) {
		vector<int> path = pathcover.at(i);
		for(int j=0;j<path.size();j++) {
			int v = path.at(j);
			this->pathsforv[v].push_back(i);
			index[v][i] = j;
		}
	}

	int last2reach[this->SGraph->getNoOfVertices()][this->pathcover.size()];

	// Init all to -1
	for(int v=0;v<this->SGraph->getNoOfVertices();v++) {
		for(int i=0;i<pathcover.size();i++)
			last2reach[v][i] = -1;
	}

	for(int v=0;v<this->SGraph->getNoOfVertices();v++) {
		std::vector<int> vpaths = pathsforv[v];
		for(int i=0;i<pathcover.size();i++) {
			// i in paths[v]
			if(find(vpaths.begin(),vpaths.end(),i)!=vpaths.end()) {
				last2reach[v][i] = index[v][i];
			}
			// i not in paths[v]
			else {
				vector<Vertex*> inneighbors = SGraph->getVertex(v)->getInNeighbors();
				for(int k=0;k<inneighbors.size();k++) {
					int u = inneighbors.at(k)->getId();
					if(last2reach[u][i] > last2reach[v][i]) {
						last2reach[v][i] = last2reach[u][i];
					}
				}
			
			}
		}
	}

	for(int v=0;v<this->SGraph->getNoOfVertices();v++) {
		for(int i=0;i<pathcover.size();i++) {
			if(last2reach[v][i] != -1) {
				forward.at(pathcover.at(i).at(last2reach[v][i])).push_back(std::make_pair(v,i));
			}
		}
	}

/*
	// How many'th vertex on its path is this vertex
	int index[this->SGraph->getNoOfVertices()][this->pathcover.size()];
	int parent[this->SGraph->getNoOfVertices()][this->pathcover.size()];

	// Initialize pathsforv vector (corresponds to paths in the article)
	// and forward
	for(int i=0;i<this->SGraph->getNoOfVertices();i++) {
		std::vector<int> tmp;
		this->pathsforv.push_back(tmp);
		std::vector<std::pair<int,int> > tmp2;
		this->forward.push_back(tmp2);
	}



	for(int i=0;i<this->pathcover.size();i++) {
		vector<int> path = pathcover.at(i);
		for(int j=0;j<path.size();j++) {
			int v = path.at(j);
			this->pathsforv[v].push_back(i);
			index[v][i] = j;
		}
	}

	for(int v=0;v<this->SGraph->getNoOfVertices();v++) {

		std::vector<int> vpaths = pathsforv[v];
		
		for(int i=0;i<pathcover.size();i++) {
			// i in paths[v]
			if(find(vpaths.begin(),vpaths.end(),i)!=vpaths.end()) {
				parent[v][i] = v;
			}
			// i not in paths[v]
			else {
				vector<Vertex*> inneighbors = SGraph->getVertex(v)->getInNeighbors();

				int maxu = -1;
				int argmaxu = -1;

				for(int k=0;k<inneighbors.size();k++) {
					int u = inneighbors.at(k)->getId();
					if(index[u][i] > maxu) {
						maxu = index[u][i];
						argmaxu = u;
					}
				}
				parent[v][i] = maxu;
				index[v][i] = maxu;
			}
		}
	}

	for(int v=0;v<this->SGraph->getNoOfVertices();v++) {
		vector<int> vpaths = pathsforv[v];
		for(int i=0;i<vpaths.size();i++) {
			this->forward.at(index[v][vpaths.at(i)]).push_back(std::make_pair(v,vpaths.at(i)));
		}

		vector<int> vpaths = pathsforv[v];
		for(int i=0;i<pathcover.size();i++) {
			 if(find(vpaths.begin(),vpaths.end(),i)!=vpaths.end()) {
				this->forward.at(index[v][i]).push_back(std::make_pair(v,i));
			}
		}
	}*/
}

// +1 every key coordinate to shift from 0-based to 1-based (otherwise the key 0 causes issues)
std::vector<Tuple*> ColinearSolver::colinearChain(std::vector<Tuple*> M) {
	RMaxQTree* ITrees = new RMaxQTree[this->pathcover.size()];
	RMaxQTree* TTrees = new RMaxQTree[this->pathcover.size()];
		
	// Get the key list from anchors: all unique M[j].d
	std::vector<int> keys;
	keys.push_back(0);
		
	for(int i=0;i<M.size();i++) {
			if(std::find(keys.begin(), keys.end(), M.at(i)->d+1) == keys.end())
				keys.push_back(M.at(i)->d+1);
	}
		
	std::sort(keys.begin(),keys.end());
	
	// Change into array for the RMaxQTree
	int keysarray[keys.size()];
	
	for(int i=0;i<keys.size();i++)
		keysarray[i] = keys.at(i);
		
	// Init the trees
	for(int i=0;i<this->pathcover.size();i++) {
		ITrees[i].fillRMaxQTree(keysarray, keys.size());
		TTrees[i].fillRMaxQTree(keysarray, keys.size());
		ITrees[i].update(0,-1,0);
		TTrees[i].update(0,-1,0);
	}

	std::vector<int>* start = new std::vector<int>[this->SGraph->getNoOfVertices()];
	std::vector<int>* end = new std::vector<int>[this->SGraph->getNoOfVertices()];
		
	for(int j=0;j<M.size();j++) {
		start[M.at(j)->PFirst].push_back(j);
		end[M.at(j)->PLast].push_back(j);
	}
	
	
	for(int v=0;v<this->SGraph->getNoOfVertices();v++) {

		// For every tuple whose path ends at this vertex, update the trees for all the paths on which this vertex lies
		for(int j=0;j<end[v].size();j++) {

			std::vector<int> pathsv = this->pathsforv[v];
			for(int i=0;i<pathsv.size();i++) {

				// If the tuple is of length 1, need to update coverage first
				if(M.at(end[v].at(j))->P.size() == 1) {
					std::pair<int,int> atuple = TTrees[pathsv.at(i)].query(0,M.at(end[v].at(j))->c+1-1);
					std::pair<int,int> btuple = ITrees[pathsv.at(i)].query(M.at(end[v].at(j))->c+1, M.at(end[v].at(j))->d+1);

					int Ca = (M[end[v].at(j)]->d-M[end[v].at(j)]->c+1) + atuple.second;
					int Cb = M[end[v].at(j)]->d + btuple.second;
					
					if(Ca > M.at(end[v].at(j))->C && Ca >= Cb) {
						M.at(end[v].at(j))->C = Ca;
						if(atuple.first != -1)
							M.at(end[v].at(j))->previous = M.at(atuple.first);
					}
					else if(Cb > M.at(end[v].at(j))->C) {
						M.at(end[v].at(j))->C = Cb;
						if(btuple.first != -1)
							M.at(end[v].at(j))->previous = M.at(btuple.first);
					}				
				}

				TTrees[pathsv.at(i)].update(M.at(end[v].at(j))->d+1, end[v].at(j),M.at(end[v].at(j))->C);
				ITrees[pathsv.at(i)].update(M.at(end[v].at(j))->d+1, end[v].at(j),M.at(end[v].at(j))->C-(M.at(end[v].at(j))->d));
			}
		}

		std::vector<std::pair<int,int> > forwardv = this->forward[v];
		
		// And then follow the forward links to update the coverages for tuples whose first node is in the forward set
		for(int f=0;f<forwardv.size();f++) {
			std::vector<int> startw = start[forwardv.at(f).first];
			
			for(int j=0;j<startw.size();j++) {

				// Skip if it's link to self and path of length 1, was processed above
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
	
	std::vector<Tuple*> solution;
	
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
	solution.push_back(currentTuple);
	
	while(currentTuple->previous != NULL) {
		currentTuple =currentTuple->previous;
		solution.push_back(currentTuple);
	}
	
	// It's backwards, so reverse
	std::vector<Tuple*> revsolution;
	
	for(int i=solution.size()-1;i>=0;i--)
		revsolution.push_back(solution.at(i));
	
	return revsolution;
	
}

// Note: Patterns have to be the same length
ColinearSolver::ColinearSolver(std::string graphFile, std::string nodesFile, int patlen) {
	SGraph = new SequenceGraph();
	SGraph->readSplicingGraphFile(graphFile, nodesFile, patlen);
}

// Solves the co-linear chaining problem (patterns are read from the file, with one pattern per row)
// Outputs to the given file with first row having the pattern, following row the number of anchors,
// and then all the anchor tuples, one per row ((patstart,patend),(path in graph))
void ColinearSolver::solve(std::string patternfile, std::string outputfile, int threshold) {

	std::string tempgraphfile = "flowgraph.tmp";

	convertSequenceGraphToFlowGraph(tempgraphfile);
	solveForPaths(tempgraphfile);

	computeForward();

	std::ifstream patin(patternfile.c_str());
	std::ofstream out(outputfile.c_str());


	// Output header that tells how many vertices there are and how many paths
	out << "##Vertices:\t" << this->SGraph->getNoOfVertices() << std::endl;
	out << "##Paths:\t" << this->pathcover.size() << std::endl;

	std::string line;

	while(getline(patin, line)) {
		std::vector<Tuple*> anchors = SGraph->findAnchors(line, threshold); 

		std::vector<Tuple*> chain = this->colinearChain(anchors);


		// Output the chain
//		for(int i=0;i<chain.size();i++)
//			std::cout << chain.at(i)->toString() << ", ";
//		std::cout << std::endl;


		for(int i=0;i<chain.size();i++)
			out << chain.at(i)->toString() << ", ";

		out << std::endl;	
	}

	// TODO possible check that the file stream is good here


	patin.close();
	out.close();

}
	



