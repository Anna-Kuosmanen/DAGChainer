/*
 *
 * SequenceGraph.cpp
 *
 * Created on: Sept 17th 2017
 *	Author: aekuosma
 *
 */

// TODO: Support for different lengths of pattern for naive anchor search? Currently the pattern length has to be initialized when graph is read from files

#include "SequenceGraph.h"
#include "Tuple.h"
#include "utils.h"
#include "Vertex.h"
#include "gcsa/gcsa.h"
#include "gcsa/lcp.h"

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

// For sorting the anchors
struct less_than_path {
	inline bool operator() (const Tuple* tuple1, const Tuple* tuple2) {
		if(tuple1->PFirst != tuple2->PFirst)
			return (tuple1->PFirst < tuple2->PFirst);
		else
			return (tuple1->PLast < tuple2->PLast);
	}
};


void SequenceGraph::addToVertices(Vertex* v) {
	(this->vertices).push_back(v);
}

// Removes duplicates from the anchor list
// TODO Should this be done while searching for anchors or is it more effective here?
void SequenceGraph::removeDuplicates(std::vector<Tuple*> &anchors) {
	std::vector<Tuple*> tmp;

	for(unsigned i=0;i<anchors.size();i++) {
		bool found = false;
		for(unsigned j=0;j<tmp.size();j++) {
			if(anchors.at(i)->equals(*(tmp.at(j)))) {
				found = true;
				break;
			}
		}
 
		if(!found)
			tmp.push_back(anchors.at(i));
		// Stop memory leak
		else
			delete anchors.at(i);
	}
	anchors = tmp;
}

// ColinearSolver can't currently handle overlaps on the paths, split anchors to remove them
// Only suffix-prefix overlaps are a problem, everything else can be ignored since they can't be part of the same chain
void SequenceGraph::removePathOverlaps(std::vector<Tuple*> &anchors) {

	std::vector<int>* pathstarts = new std::vector<int>[this->getNoOfVertices()];
	std::vector<int>* pathends = new std::vector<int>[this->getNoOfVertices()];

	// For every vertex, mark which paths start and end here
	for(unsigned i=0;i<anchors.size();i++) {
		pathstarts[anchors.at(i)->PFirst].push_back(i);
		pathends[anchors.at(i)->PLast].push_back(i);
	}

	std::vector<Tuple*> fixed_anchors;
	
	// For every path, check if some other path starts in it
	for(unsigned i=0;i<anchors.size();i++) {
		Tuple* tup = anchors.at(i);
		bool modified = false;
		int last_start = tup->c;
		int last_path_start = tup->PFirst;

		std::vector<int> new_path;
		new_path.push_back(tup->PFirst);
		// Containments are ok, so if they start on same node, it doesn't matter
		for(unsigned j=1;j<tup->P.size();j++) {
			// Something starts here, cut right before it
			if(pathstarts[tup->P.at(j)].size() > 0) {

				std::vector<int> starts = pathstarts[tup->P.at(j)];
				modified = true;

				fixed_anchors.push_back(new Tuple(new_path, last_start, last_start+new_path.size()-1));
				last_start = last_start+new_path.size();
				new_path.clear();
			}
			new_path.push_back(tup->P.at(j));
		}
		// And then add the remaining part
		if(modified && last_start <= tup->d) {
			fixed_anchors.push_back(new Tuple(new_path, last_start, last_start+new_path.size()-1));
		}


		// Clean up if this tuple isn't included after fix
		if(modified)
			delete anchors.at(i);
		else
			fixed_anchors.push_back(anchors.at(i));
	}


	anchors = fixed_anchors;


	delete[] pathstarts;
	delete[] pathends;
}

SequenceGraph::SequenceGraph() {
}

// Since we use Vertex pointers, need to delete them explicitly here
SequenceGraph::~SequenceGraph() {
	for(int i=0;i<this->getNoOfVertices();i++)
		delete this->vertices.at(i);

	vertices.clear();
}

int SequenceGraph::getNoOfVertices() {
	return this->vertices.size();
}

// Returns i'th vertex (or NULL if invalid)
Vertex* SequenceGraph::getVertex(int i) {
	if(i < int(this->vertices.size()))
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

	for(unsigned j=0;j<neighbors.size();j++) {
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

void SequenceGraph::clearDArray(int patlen) {
	for(int v=0;v<this->getNoOfVertices();v++) {
		for(int i=0;i<patlen;i++) {
			this->getVertex(v)->updateDValue(i, 0);
		}
	}
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
		// TODO Is this correct? Why report range of length 1?
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

void SequenceGraph::backtrackPath(std::vector<int> &path, std::string seq, int curpos, Vertex* curvertex) {

	// If path is of same length as seq, some other branch found the whole thing, return
	if(path.size() == seq.size())
		return;
	
	// No match, return
	if(curvertex->getLabel() != seq.at(curpos))
		return;

	// Otherwise if the label matched, add your id to path
	path.push_back(curvertex->getId());

	// If whole seq hasn't been used, call recursively for your outneighbors
	if(curpos < seq.size()-1) {

		std::vector<Vertex*> neighbors = curvertex->getOutNeighbors();

		for(unsigned i=0;i<neighbors.size();i++)
			this->backtrackPath(path, seq, curpos+1, neighbors.at(i));

	

		// If after recursing on all outneighbors the sequence hasn't been used up, this couldn't have been it, pop
		if(path.size() != seq.size()) {
			
			while(path.at(path.size()-1) != curvertex->getId()) {
				path.pop_back();
			}
			// And self
			path.pop_back();

		}
	}

}

// Going backwards from end
void SequenceGraph::backtrackReversePath(std::vector<int> &path, std::string seq, int curpos, Vertex* curvertex) {
	// If path is of same length as seq, some other branch found the whole thing, return
	if(path.size() == seq.size())
		return;
	
	// No match, return
	if(curvertex->getLabel() != seq.at(curpos))
		return;

	// Otherwise if the label matched, add your id to path
	path.push_back(curvertex->getId());

	// If whole seq hasn't been used, call recursively for your inneighbors
	if(curpos > 0) {

		std::vector<Vertex*> neighbors = curvertex->getInNeighbors();
		for(unsigned i=0;i<neighbors.size();i++)
			this->backtrackReversePath(path, seq, curpos-1, neighbors.at(i));

		// If after recursing on all outneighbors the sequence hasn't been used up, this couldn't have been it, pop
		if(path.size() != seq.size()) {
			
			while(path.at(path.size()-1) != curvertex->getId()) {
				path.pop_back();
			}
			// And self
			path.pop_back();

		}
	}

}



// Need to do this here, since Node::decode only gives the first node
// Need to backtrack with the sequence on the sequencegraph
// TODO This only works right if GFA is created with "concat" option! Fix.
void SequenceGraph::convertMEMToTuple(std::vector<Tuple*> &tuples, std::vector<Tuple*> &revtuples, MaximalExactMatch mem, std::string seq) {

	// Each node is a startpoint of a path
	// For each node, track through the sequence graph with the pattern
	for(auto& node : mem.nodes) {
		// This is of form id of node : (-) offset (where - means reverse complement)
		std::string nodecode = gcsa::Node::decode(node);
		std::vector<std::string> parts = split(nodecode, ':');
		bool rc = false;

		if(parts.at(1)[0] == '-') {
			rc = true;
			parts.at(1) = parts.at(1).substr(1,std::string::npos);
		}

		std::vector<int> path;


		// Reversecomplements are a little more complicated
		if(rc) {

			// The "first" (last) node is offset from the end of the exon
			// -1 on the splicing nodes because had to +1 all coordinates for vg graph
			// Sequence needs to be reverse complemented, and start from the end
			int nodenumber = this->splicingNodesToSequenceNodes.at(std::atoi(parts.at(0).c_str())-1).second-(std::atoi(parts.at(1).c_str()));
			Vertex* curvertex = this->getVertex(nodenumber);
			this->backtrackReversePath(path, reverseComplement(seq), int(seq.size()-1), curvertex);
		}
		else {
			// The node number of sequence graph is the first node of that exon + offset
			int nodenumber = this->splicingNodesToSequenceNodes.at(std::atoi(parts.at(0).c_str())-1).first+std::atoi(parts.at(1).c_str());

	
			Vertex* curvertex = this->getVertex(nodenumber);

			this->backtrackPath(path, seq, 0, curvertex);
		}
	
		if(rc)
			revtuples.push_back(new Tuple(path, mem.begin, mem.end));
		else
			tuples.push_back(new Tuple(path, mem.begin, mem.end));

	}

}


// Find all the anchors between a pattern and the sequence graph
// NOTE: This is an old naive method, use findAnchorsGSCA2 for speed
void SequenceGraph::findAnchorsNaive(std::vector<Tuple*> &anchors, std::string pattern, int threshold) {

	// For every vertex
	for(int ver=0;ver<this->getNoOfVertices();ver++) {
		// For every length of prefix
		for(unsigned i=0;i<=pattern.size();i++) {
			bool retval = this->fillDArray(this->vertices.at(ver),i, pattern);
			if(!retval) {
				// Label of v didn't match the pattern for this i
				// Report on every in-neighbor of v for index i-1, if D value is over 0
				std::vector<Vertex*> inneighbors = this->vertices.at(ver)->getInNeighbors();

				for(unsigned j=0;j<inneighbors.size();j++) {
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
	this->clearDArray(pattern.size());

}

void SequenceGraph::findAnchorsGCSA2(std::vector<Tuple*> &anchors, std::vector<Tuple*> &revanchors, std::string gcsa_name, std::string lcp_name, std::string pattern, int threshold) {
	gcsa::GCSA gcsa;
	sdsl::load_from_file(gcsa, gcsa_name);

	gcsa::LCPArray lcp;
	sdsl::load_from_file(lcp, lcp_name);

	std::vector<MaximalExactMatch> mems;

	// Empty pattern, don't return anything (it'd match everywhere)
	if(pattern.size() == 0) {
		std::vector<Tuple*> tuples;
		return;
	}
	// Starting from the end
	int current_pos = pattern.size()-1;

	auto full_range = gcsa::range_type(0, gcsa.size()-1);
	MaximalExactMatch match = MaximalExactMatch(current_pos, current_pos, full_range);


	gcsa::range_type last_range = match.range;
	while(current_pos >= 0) {
		// Need to remember this, if currently processed char doesn't match
		last_range = match.range;

		match.range = gcsa.LF(match.range, gcsa.alpha.char2comp[pattern.at(current_pos)]);

		// No match or exceeded the order of index, report and go to parent to chop from end
		if(gcsa::Range::empty(match.range) || match.end-current_pos > gcsa.order()) {
			match.begin = current_pos+1;
			match.range = last_range;
			// Report if length exceeds the threshold
			// The first check blocks placeholder if the threshold=0 for some reason
			if(match.end >= match.begin && match.end-match.begin+1 >= threshold) {
				mems.push_back(match);

			}
			// Note: this is following the example of VG's mapper.cpp, not sure what these step sizes mean
			size_t last_mem_length = match.end - match.begin+1;
			gcsa::STNode parent = lcp.parent(last_range);
			size_t step_size = last_mem_length - parent.lcp();
			match.end = match.end-step_size;
			match.range = parent.range();
		}
		// Matching
		else {
			match.begin = current_pos;
			current_pos--;
		}

	}

	// Checking if there's a MEM (last char matched, so it didn't go to above loop's if)
	if(match.end-match.begin > threshold) {
		mems.push_back(match);
	}

	// Convert to tuples
	for(unsigned i=0;i<mems.size();i++) {
		gcsa.locate(mems.at(i).range, mems.at(i).nodes);

		this->convertMEMToTuple(anchors, revanchors, mems.at(i), pattern.substr(mems.at(i).begin, mems.at(i).end-mems.at(i).begin+1));
	}

	// Reverse anchors give the positions in the original pattern that match the reverse complement of the graph
	// -> have to fix the positions in the pattern to match the reverse complement of the pattern
	// (otherwise we'd have to look for a colinear chain that goes in descending order in the pattern in ColinearSolver)
	for(unsigned i=0;i<revanchors.size();i++) {
		Tuple* anchor = revanchors.at(i);
		int newc = pattern.size()-anchor->d-1;
		int newd = pattern.size()-anchor->c-1;
		anchor->c = newc;
		anchor->d = newd;

	}

	sort(anchors.begin(), anchors.end(), less_than_path());
	sort(revanchors.begin(), revanchors.end(), less_than_path());

	// Remove path overlaps (ColinearSolver implementation currently doesn't allow for them, would need to add two-dimensional range queries)
	this->removePathOverlaps(anchors);
	this->removePathOverlaps(revanchors);

}


// Graph file has first row the number of nodes (exons), then the edges (i+1'th row has the outgoing edges from i'th node)
// Nodes file has pairs of node number (first row) and the corresponding sequence (second row) for building the sequence graph (they need to be in order, labels are not checked currently!)
// TODO Check the labels, then the nodes file wouldn't need to be in order
// Patlen is the length of the reads (needed to initialize the arrays in vertex objects)
// Returns true if the reading succeeded (files existed)
bool SequenceGraph::readSplicingGraphFile(std::string graphfile, std::string nodesfile, int patlen) {
	std::ifstream graphIn(graphfile.c_str());
	std::ifstream nodesIn(nodesfile.c_str());

	if(!graphIn.good() || !nodesIn.good()) {
		graphIn.close();
		nodesIn.close();
		return false;

	}

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
			for(unsigned j=0;j<content.size();j++) {
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
			for(unsigned j=0;j<parts.size();j++) {
				iparts.push_back(std::atoi(parts.at(j).c_str()));
			}

			for(unsigned j=0;j<iparts.size();j++) {
				this->vertices.at(lastvertices.at(i))->addEdgeTo(this->vertices.at(firstvertices.at(iparts.at(j))));
			}
		}
	}

	// Go over the firstvertices and lastvertices to save the mapping from splicing graph nodes to sequence graph nodes
	for(unsigned i=0;i<firstvertices.size();i++) {
		splicingNodesToSequenceNodes.push_back(std::make_pair(firstvertices.at(i),lastvertices.at(i)));

	}


	graphIn.close();
	nodesIn.close();
	return true;
}

// TODO Splicing graphs are by default all in forward direction
// If/when add variation graph, need to consider inverses (+/- or -/+)
// as well as overlaps (last column)
// TODO Sanity check, that all sequence nodes match to some exon node
void SequenceGraph::outputGFA(std::string gfaout, bool concat) {

	std::ofstream out(gfaout.c_str());

	out << "H\tVN:Z:1.0" << std::endl;

	// VG doesn't accept node index 0, have to +1 all

	if(!concat) {
		for(int i=0;i<this->getNoOfVertices();i++) {

			// The vertex as segment
			out << "S\t" << i+1 << "\t" << this->getVertex(i)->getLabel() << std::endl;

			// All the outneighbors as links
			std::vector<Vertex*> outneighbors = this->getVertex(i)->getOutNeighbors();

			for(unsigned j=0; j<outneighbors.size();j++) {

				out << "L\t" << i+1 << "\t+\t" << outneighbors.at(j)->getId()+1 << "\t+\t0" << std::endl;

			}
		}
	}
	else {

		std::stringstream sstm;

		// The original node coordinates are saved in splicingNodesToSequenceNodes
		for(unsigned i=0;i<splicingNodesToSequenceNodes.size();i++) {

			// Concatenate the labels to form a segment and print
			for(int j=splicingNodesToSequenceNodes.at(i).first;j<=splicingNodesToSequenceNodes.at(i).second;j++) {
				sstm << this->getVertex(j)->getLabel();

			}
			out << "S\t" << i+1 << "\t" << sstm.str() << std::endl;
			sstm.str("");


			// Check the outneighbors of the last node
			std::vector<Vertex*> outneighbors = this->getVertex(this->splicingNodesToSequenceNodes.at(i).second)->getOutNeighbors();

			for(unsigned n=0;n<outneighbors.size();n++) {
				int neighborID = outneighbors.at(n)->getId();

				// Find for which exon this is the first character
				for(unsigned s=i+1;s<splicingNodesToSequenceNodes.size();s++) {
					if(splicingNodesToSequenceNodes.at(s).first == neighborID) {
						out << "L\t" << i+1 << "\t+\t" << s+1 << "\t+\t0" << std::endl;
						break;
					}
				}
			}
		}
		
	}

	out.close();

}

// Helper function for convertChainToExons
// Colinear chain might not have any anchors on short exons, so here we check that the path is valid, and add nodes if necessary
// If adding more than one exon between i'th and i+1'th exon of the path would be required, or adding more than "stringency" parameter exons total, scrap the whole thing, it's too unreliable
void SequenceGraph::fixPath(std::vector<int> &exons, int stringency) {

	std::vector<int> newexons;
	int addedCount = 0;
	// Flag as true if could not find an exon between i'th and i+1'th that would fix the path
	bool toohard = false;

	for(unsigned i=0;i<exons.size()-1;i++) {

		bool arcok = false;

		std::vector<Vertex*> outneighbors = this->getVertex(splicingNodesToSequenceNodes.at(exons.at(i)).second)->getOutNeighbors();

		for(unsigned j=0;j<outneighbors.size();j++) {
			if(outneighbors.at(j)->getId() == splicingNodesToSequenceNodes.at(exons.at(i+1)).first) {
				arcok = true;
				break;
			}
		}

		newexons.push_back(exons.at(i));

		// There was no edge between i'th and i+1'th exon, seek for addition
		// Ties are broken in favor of the exon with smaller index
		if(!arcok) {
			int startExon = exons.at(i);
			int endExon = exons.at(i+1);

			std::vector<Vertex*> startNeighbors = this->getVertex(splicingNodesToSequenceNodes.at(startExon).second)->getOutNeighbors();
			std::vector<Vertex*> endNeighbors = this->getVertex(splicingNodesToSequenceNodes.at(endExon).first)->getInNeighbors();

			bool found = false;

			for(unsigned k=startExon+1;k<endExon;k++) {
				bool startNeighborFound = false;
				bool endNeighborFound = false;

				for(unsigned n=0;n<startNeighbors.size();n++)
					if(startNeighbors.at(n)->getId() == splicingNodesToSequenceNodes.at(k).first)
						startNeighborFound = true;

				for(unsigned n=0;n<endNeighbors.size();n++)
					if(endNeighbors.at(n)->getId() == splicingNodesToSequenceNodes.at(k).second)
						endNeighborFound = true;

				if(startNeighborFound && endNeighborFound) {
					found = true;
					newexons.push_back(k);
					break;
				}

			}

			if(!found) {
				toohard = true;
				break;
			}
			else
				addedCount++;

		}



	}

	newexons.push_back(exons.at(exons.size()-1));

	if(addedCount > stringency || toohard)
		newexons.clear();

	exons = newexons;

}

std::vector<int> SequenceGraph::convertChainToExons(ColinearChain chain, int stringency) {

	std::vector<Tuple*> tuples = chain.chain;
	std::vector<int> exons;

	// Find the exon where this path starts
	int current_exon = -1;

	for(unsigned j=0;j<splicingNodesToSequenceNodes.size();j++) {
		std::pair<int,int> exon = splicingNodesToSequenceNodes.at(j);

		if(tuples.at(0)->PFirst >= exon.first && tuples.at(0)->PFirst <= exon.second) {
			current_exon = j;
			exons.push_back(current_exon);
			break;
		}
	}



	for(unsigned i=1;i<tuples.size();i++) {
		
		std::vector<int> path = tuples.at(i)->P;

		for(unsigned j=1;j<path.size();j++) {

			if(path.at(j) >= splicingNodesToSequenceNodes.at(current_exon).first && path.at(j) <= splicingNodesToSequenceNodes.at(current_exon).second) {
				continue;
			}

			// If go over the end of current exon, find the next exon
			else {
				for(unsigned k=current_exon;k<splicingNodesToSequenceNodes.size();k++) {
					std::pair<int,int> exon = splicingNodesToSequenceNodes.at(k);

					if(path.at(j) >= exon.first && path.at(j) <= exon.second) {
						current_exon = k;
						exons.push_back(current_exon);
						break;
					}
				}
			}
		}
	}

	this->fixPath(exons, stringency);

	return exons;

}

