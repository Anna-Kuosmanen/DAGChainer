/*
 * MPCGraph.h
 *
 *  Created on: Feb 11, 2014
 *      Author: ahmedsobih and aekuosma
 */

#ifndef MPCGRAPH_H_
#define MPCGRAPH_H_
#include "MPCHeaders.h"
#include "MPCNode.h"
#include "MPCArc.h"
#include "MPCLogger.h"
#include "MPCUtil.h"
#include <tuple>

typedef tuple<int,int,int> overlaptuple;


bool mycompare (const overlaptuple &lhs, const overlaptuple &rhs){
	return get<2>(lhs) < get<2>(rhs);
}


class MPCGraph {
private:
	map<int,MPCNode> nodeMap;
	map<int, MPCArc> arcMap;
	vector< vector<int> > subpathConstraintVector;
	set<int> startSolNodeSet;
	set<int> endSolNodeSet;
	string solution;
	double maxConverage;
	//correspondence between nodes in flowNetwork and nodes in g
	map<int, MPCNode> flowNodeToMPCNode;
	map<int, MPCArc> flowArcToMPCArc;

	void createSolPath(ListDigraph& flowNetwork, ListDigraph::ArcMap<int64_t>& flowMap, ListDigraph::Node flNode, string& solution, int& totalFlow){
		MPCNode& mpcSourceNode=flowNodeToMPCNode.find(flowNetwork.id(flNode))->second;
		if(mpcSourceNode.getId()==MPCNode::SOURCE_NODE_ID){
//			solution+="-1 ";
		}
		else if(mpcSourceNode.getId()==MPCNode::SINK_NODE_ID){
			solution+="\n";
		}
		for (ListDigraph::OutArcIt arc(flowNetwork,flNode); arc != INVALID; ++arc){
			if (flowMap[arc] != 0 ){
				ListDigraph::Node targetNode=flowNetwork.target(arc);
				MPCNode& mpcTargetNode=flowNodeToMPCNode.find(flowNetwork.id(targetNode))->second;
				MPCArc& mpcArc=flowArcToMPCArc.find(flowNetwork.id(arc))->second;
				if(mpcSourceNode.getId()==MPCNode::SOURCE_NODE_ID && mpcTargetNode.getId()==MPCNode::SINK_NODE_ID){
					totalFlow-=flowMap[arc];
					flowMap[arc]=0;
					continue;
				}else{
					flowMap[arc]--;
					totalFlow--;
				}
				if(mpcArc.getSubpathId()!=-1){
					vector<int> subpath=subpathConstraintVector[mpcArc.getSubpathId()];
					for(unsigned i=1; i<subpath.size();i++){
						solution+=MPCUtil::getStringValue(subpath[i])+" ";
					}
				}
				else{
					if(mpcTargetNode.getId()>MPCNode::SINK_NODE_ID)
					solution+=MPCUtil::getStringValue(mpcTargetNode.getId())+" ";
				}
				createSolPath(flowNetwork, flowMap, targetNode, solution, totalFlow);
				break;
			}
		}
	}

	// Note: the constraints MUST be sorted
	void setSubpathConstaints(vector<string>& subPathVector, vector<double>& subPathCoverages){

		if(subPathVector.size() == 0)
			return;
		set<int> pathToRemoveIdx;
		string pathA, pathB;
		// Check if one path is completely included in the other path
		for(unsigned i=0;i<subPathVector.size()-1;i++){
			pathA=subPathVector[i];
			for(unsigned j=i+1;j<subPathVector.size();j++){
				pathB=subPathVector[j];
				if(MPCUtil::isIncluded(pathA, pathB)){
					if(pathA.size()>=pathB.size())
						subPathVector[j]="";
					else
						subPathVector[i]="";
				}
			}
		}

		vector<string> temp_vector;

		for(unsigned i=0;i<subPathVector.size();i++){
//			cout << "i=" << i << ", subpath=" << subPathVector.at(i) << endl;
			if(subPathVector.at(i)!="")
				temp_vector.push_back(subPathVector.at(i));

		}

		subPathVector.clear();

//		cout << "SubpathVector:" << endl;
		for(unsigned i=0;i<temp_vector.size();i++) {
//			cout << temp_vector.at(i) << endl << endl;;
			subPathVector.push_back(temp_vector.at(i));
		}

		temp_vector.clear();

		// Create a graph of the constraints and let flow decide which paths should be merged
		vector<vector<int> > arcs;

		// Put the start and end nodes to vectors that don't have to keep splitting the paths
		vector<int> pathStartNodes;
		vector<int> pathEndNodes;

		for(unsigned i=0;i<subPathVector.size();i++) {
			vector<string> parts;
			split(parts, subPathVector.at(i), is_any_of(MPCUtil::SPACE_SEPARATOR));
			pathStartNodes.push_back(MPCUtil::getIntValue(parts[0]));
			pathEndNodes.push_back(MPCUtil::getIntValue(parts[parts.size()-1]));
		}

		for(unsigned i=0;i<subPathVector.size();i++) {
	
			vector<int> nodeArcs;

			// Check overlaps
			for(unsigned j=i+1;j<subPathVector.size();j++) {
				if(MPCUtil::mergePath(subPathVector.at(i),subPathVector.at(j)) !="") 
					nodeArcs.push_back(j);

			}	

/*			// Do traversal of the graph to check which path start nodes are reachable from this path's end node, these are the transitory arcs
			vector<int> stack;
			stack.push_back(pathEndNodes.at(i));

			bool visitedNodeIds[countNodes()];
			for(int i=0;i<countNodes();i++)
				visitedNodeIds[i]=false;

			double cost = 0.0;

			while(stack.size() > 0) {
				int node = stack.back();
				stack.pop_back();
				if(!visitedNodeIds[node]) {
					visitedNodeIds[node]=true;
					vector<MPCArc> outArcs = getOutArcs(getNode(getCopyNodeId(node)));

					for(int arcNo=0;arcNo<outArcs.size();arcNo++) {
						MPCNode neighborNode = getDestination(outArcs.at(arcNo));
						stack.push_back(neighborNode.getId());
					}

					

					for(int j=i+1;j<subPathVector.size();j++) {
						if(pathStartNodes.at(j) == -1)
							continue;
						if(pathStartNodes.at(j) == node) {
							nodeArcs.push_back(j);
						}

					}
				}

			}
*/
			arcs.push_back(nodeArcs);
		}
		vector<int> sourceConstraints;
		vector<int> sinkConstraints;

		int no_of_constraints = subPathVector.size();
		bool* source = new bool[no_of_constraints];
		bool* sink = new bool[no_of_constraints];

		for(int i=0;i<no_of_constraints;i++){
			source[i] = true;
			sink[i] = true;
		}

		for(unsigned i=0;i<arcs.size();i++) {
			// If it has outneighbors, it's not sink
			if(arcs.at(i).size() > 0)
				sink[i] = false;

		}
		// If someone points to it (inneighbor), it's not a source
		for(unsigned i=0;i<arcs.size();i++) {
			for(unsigned j=0;j<arcs.at(i).size();j++)
				source[arcs.at(i).at(j)] = false;
		}

		// But, source constraints are the subpath constraints whose first node is a source, and sinks the ones whose last node is a sink (trumps the arcs)
		set<int> sinkNodes = getEndSolNodeSet();
		set<int> sourceNodes = getStartSolNodeSet();
		for(set<int>::iterator sit=sinkNodes.begin();sit!=sinkNodes.end();sit++) {
			for(unsigned i=0;i<pathEndNodes.size();i++) {
				if(*sit == getCopyNodeId(pathEndNodes.at(i))) {
					sink[i] = true;;
				}
			}
		}
		for(set<int>::iterator sit=sourceNodes.begin();sit!=sourceNodes.end();sit++) {
			for(unsigned i=0;i<pathStartNodes.size();i++)
				if(*sit == pathStartNodes.at(i))
					source[i] = true;
		}


		for(int i=0;i<no_of_constraints;i++) {
			if(source[i])
				sourceConstraints.push_back(i);
			if(sink[i])
				sinkConstraints.push_back(i);
		}

		delete[] source;
		delete[] sink;


		MPCGraph constraintGraph = createMPCGraphFromConstraints(subPathVector, subPathCoverages, arcs, sourceConstraints, sinkConstraints);
		string constraintSolution = constraintGraph.solve();
		

//		cerr << "Constraint solution:" << endl;
//		cerr << constraintSolution << endl;


		vector<string> paths;
		split(paths, constraintSolution, is_any_of("\n"));

		// paths.size()-1 because there's an extra line feed at end
		for(unsigned i=0;i<paths.size()-1;i++) {
			vector<int> subPath;
			vector<string> path;
			split(path, paths.at(i), is_any_of(MPCUtil::SPACE_SEPARATOR));

			// First node is a dummy for forcing min number of paths,
			vector<string> new_path_vector;
	                // Dummy at start and extra space at end, so if the size is 3, there's only one node 
			if(path.size() == 3) {
				string temp_path = subPathVector.at(MPCUtil::getIntValue(path.at(1)));
				split(new_path_vector, temp_path, is_any_of(MPCUtil::SPACE_SEPARATOR));
				for(unsigned j=0;j<new_path_vector.size();j++)
					subPath.push_back(MPCUtil::getIntValue(new_path_vector.at(j)));
				this->subpathConstraintVector.push_back(subPath);
	                    continue;
	                }
			string new_path = subPathVector.at(MPCUtil::getIntValue(path.at(1)));
			// Again there's an extra space at end
			for(unsigned j=2;j<path.size()-1;j++) {
				int constraint = MPCUtil::getIntValue(path.at(j));
				string temp_path = MPCUtil::mergePath(new_path,subPathVector.at(constraint));
				if(temp_path != "")
					new_path = temp_path;
				// The else below is remnant from trying transitive closure of the graph, could be deleted
				else {
					split(new_path_vector, new_path, is_any_of(MPCUtil::SPACE_SEPARATOR));
					for(unsigned j=0;j<new_path_vector.size();j++)
						subPath.push_back(MPCUtil::getIntValue(new_path_vector.at(j)));
					this->subpathConstraintVector.push_back(subPath);
					new_path = subPathVector.at(constraint);
				}
			}
			split(new_path_vector, new_path, is_any_of(MPCUtil::SPACE_SEPARATOR));
			for(unsigned j=0;j<new_path_vector.size();j++)
				subPath.push_back(MPCUtil::getIntValue(new_path_vector.at(j)));
			this->subpathConstraintVector.push_back(subPath);
		}

		// i, j, len
/*		vector<overlaptuple> overlaps;

		// Construct compressed suffix array from the concatenation of the subpath constraints
		string concat_paths = "";
		char path_separator = (char)255;
		// Offset to make chars not be special characters
		int char_offset = 33;

		for(int i=0;i<subPathVector.size();i++){
			if(subPathVector.at(i) == "") {
				concat_paths += path_separator;
				continue;
			}
			vector<string> splitted;
			split(splitted,subPathVector.at(i),is_any_of(MPCUtil::SPACE_SEPARATOR));
			for(int j=0;j<splitted.size();j++) {
				if(MPCUtil::getIntValue(splitted.at(j))+char_offset > 254) {
					cerr << "Error! Alphabet is too large for char! Change to int alphabet (somehow). Exiting." << endl;
					exit(1);
				}
				concat_paths += (char)(MPCUtil::getIntValue(splitted.at(j))+char_offset);
			}
			concat_paths += path_separator;
		}

                //cout << "concatenated: " << endl <<  concat_paths << endl;
		

		sdsl::csa_sada<> fm_index;
		sdsl::construct_im(fm_index, concat_paths.c_str(),1);


		// Create bitvector for the concatenation, that know how many'th subpath overlaps
		sdsl::bit_vector b(concat_paths.length(),0);

		for (size_t i = 0; i < concat_paths.length(); i++) {
			if(concat_paths.at(i) == path_separator)
				b[i] = 1;
		}

		sdsl::rank_support_v<1> b_rank(&b);



		// Find the path whose prefix has the longest overlap for every path's suffix
		for(int i=0;i<subPathVector.size();i++) {

			if(subPathVector.at(i)=="")
				continue;
		//cout << "Subpathvector: " << subPathVector.at(i) << " (i=" << i << ")" << endl;

			vector<string> query;

			split(query,subPathVector.at(i),is_any_of(MPCUtil::SPACE_SEPARATOR));

			size_t l = 0;
			size_t r = fm_index.size()-1;
			size_t l2;
			size_t r2;

			int overlap_length = 0;

			for(int j=query.size()-1;j>=0;j--) {

				char query_node = (char)(MPCUtil::getIntValue(query.at(j))+char_offset);
				if(sdsl::backward_search(fm_index,l,r,query_node,l2,r2)>1) {
					l=l2;
					r=r2;
					overlap_length++;
				}
				else {
					break;
				}
			}

			// Search for path separator that would signal for prefix
			if(sdsl::backward_search(fm_index,l,r,path_separator,l2,r2)>0 && overlap_length > 0) {
				// Remember to +1 rank, because it doesn't count the current character ("1s _before_ this position")
				for(int k=l2;k<=r2;k++) {
					int overlap_index = b_rank(fm_index[k]+1);
					if(overlap_index != i)
						overlaps.push_back(make_tuple(i,overlap_index,overlap_length));
				}
			}

		}

		// Sort the tuples based on the length of the overlap
		sort(overlaps.begin(),overlaps.end(),mycompare);



		for(int i=0;i<overlaps.size();i++)
			cerr << "Overlap " << get<0>(overlaps.at(i)) << ", " << get<1>(overlaps.at(i))<< ", " << get<2>(overlaps.at(i)) << endl;


		// Create two boolean arrays to tell if each path has been used as suffix of prefix in merging
		bool *mergedAsPrefix = new bool[subPathVector.size()];
		bool *mergedAsSuffix = new bool[subPathVector.size()];
		for(int i=0;i<subPathVector.size();i++) {
			mergedAsPrefix[i] = false;
			mergedAsSuffix[i] = false;
		}

		// Represent paths as double-linked lists (see Rizzi, Tomescu and Mäkinen. RECOMB-Seq 2014)
		doubleLinkedList<int>* paths_as_lists[subPathVector.size()];

		for(int i=0;i<subPathVector.size();i++) {
			vector<string> splitted;
			split(splitted, subPathVector.at(i), is_any_of(MPCUtil::SPACE_SEPARATOR));
			vector<int> splitted_ints;
			for(int j=0;j<splitted.size();j++)
				splitted_ints.push_back(MPCUtil::getIntValue(splitted.at(j)));
			paths_as_lists[i] = new doubleLinkedList<int>(splitted_ints, splitted_ints.size());
		}

		// Start merging
		while(overlaps.size() > 0) {
			overlaptuple tup = overlaps.back();
			overlaps.pop_back();
			int first_index = get<0>(tup);
			int second_index = get<1>(tup);	
			int overlap_length = get<2>(tup);
cerr << first_index << ", " << second_index << ", " << overlap_length << endl;

			if(MPCUtil::isIncluded(subPathVector.at(first_index), subPathVector.at(second_index)))
				continue;

			vector<string> first_parts;
			split(first_parts, subPathVector.at(first_index), is_any_of(MPCUtil::SPACE_SEPARATOR));
			vector<string> second_parts;
			split(second_parts, subPathVector.at(second_index), is_any_of(MPCUtil::SPACE_SEPARATOR));

			int first_node_non_overlap = MPCUtil::getIntValue(first_parts.at(first_parts.size()-1-overlap_length));
			int second_node_non_overlap = MPCUtil::getIntValue(second_parts.at(overlap_length));

			if(!(mergedAsPrefix[first_index]) && !(mergedAsSuffix[second_index])) {

				mergedAsPrefix[first_index] = true;
				mergedAsSuffix[second_index] = true;
				doubleLinkedList<int>::mergeLists(paths_as_lists[first_index],paths_as_lists[second_index],overlap_length);
			}
		}
		// Get the merged paths and push them to subpathConstraintVector
		for(int i=0;i<subPathVector.size();i++) {
			vector<int> subPath;
			doubleLinkedList<int>* path = paths_as_lists[i];

			// Either list is empty or first node isn't start node, so don't count this path
			if(path->getHead() == NULL)
				continue;
			if(path->getHead()->prev != NULL)
				continue;
			doubleLinkedList<int>::node* p = path->getHead();
			while(p != NULL) {
				subPath.push_back(p->data);
				p=(p->next);
			}
			this->subpathConstraintVector.push_back(subPath);
			delete path;
		}

		delete[] mergedAsPrefix;
		delete[] mergedAsSuffix;

		//for(int i=0;i<subpathConstraintVector.size();i++) {
			//cout << "constraints vector:" << endl;
			//vector<int> subpath = subpathConstraintVector.at(i);
			//for(int j=0;j<subpath.size();j++)
				//cout << subpath.at(j) << " ";
			//cout << endl;
		//}
*/
/*		vector<string> pathNodes;
		bool notMerged=true;
		string mergedPath;
		while(notMerged){
			notMerged=false;
			for(int i=0;i<subPathVector.size()-1;i++){
				pathA=subPathVector[i];
				if(pathA=="")
					continue;
				for(int j=i+1;j<subPathVector.size();j++){
					pathB=subPathVector[j];
					if(pathB=="")
						continue;
					mergedPath=MPCUtil::mergePath(pathA,pathB);
					if(mergedPath!=""){
						subPathVector[i]=pathA=mergedPath;
						subPathVector[j]="";
						notMerged=true;
					}
				}
			}
		}
		for(int i=0;i<subPathVector.size();i++){
			vector<int> subPath;
			if(subPathVector[i]=="")
				continue;
			split(pathNodes, subPathVector[i], is_any_of(MPCUtil::SPACE_SEPARATOR));
			for(int j=0;j<pathNodes.size();j++){
				subPath.push_back(MPCUtil::getIntValue(pathNodes[j]));
			}
			this->subpathConstraintVector.push_back(subPath);
		}*/
	}
	MPCGraph(){
		maxConverage=0;
	}
	void createMPCGraph(string fileName){
		ifstream graphFile(fileName.c_str());
		if (!graphFile.is_open()) { // check for successful opening
			MPCLogger::log("File "+fileName+" does not exist!");
			return ;
		}
		vector<string> splittedLine;
		string line;
		/*********Create Nodes*********/
		readLine(graphFile,line); // read Number of nodes

		int numOfNodes=MPCUtil::getIntValue(line);
		for (int i=0;i<numOfNodes;i++){
			addNode(i); // create Nodes
			addCopyNode(i);
		}
		MPCLogger::log(MPCUtil::getStringValue(numOfNodes)+" nodes were added.");
		/*********Create Arcs*********/
		int arcId=0;
		for(int i=0;i<numOfNodes;i++){
			getline(graphFile,line);
			if(line=="")
				continue;
			split(splittedLine, line, is_any_of(MPCUtil::SPACE_SEPARATOR));
			for(unsigned j=0;j<splittedLine.size();j++){
				addArc(arcId, getCopyNodeId(i), MPCUtil::getIntValue(splittedLine[j])); // Create arc
				arcId++;
			}
		}
		for(int i=0;i<numOfNodes;i++){
			addArc(arcId, i, getCopyNodeId(i)); // Create arc
			arcId++;
		}
		MPCLogger::log(MPCUtil::getStringValue(arcId)+" arcs were added.");
		/*********Set Nodes Coverage*********/
		readLine(graphFile,line); // read Node coverage line
		split(splittedLine, line, is_any_of(MPCUtil::SPACE_SEPARATOR));
		double coverage;
		for(unsigned j=0;j<splittedLine.size();j++){
			coverage=MPCUtil::getDoubleValue(splittedLine[j]);
			getNode(j).setCoverage(coverage);
			getNode(getCopyNodeId(j)).setCoverage(coverage);
			getArc(j,getCopyNodeId(j)).setCoverage(coverage);
			if(maxConverage<coverage)
				maxConverage=coverage;
		}

		/*********Set Arcs Coverage*********/
		arcId=0;
		for(int i=0;i<numOfNodes;i++){
			readLine(graphFile,line); // Read arcs coverage line
			if(line=="")
				continue;
			split(splittedLine, line, is_any_of(MPCUtil::SPACE_SEPARATOR));
			for(unsigned j=0;j<splittedLine.size();j++){
				coverage=MPCUtil::getDoubleValue(splittedLine[j]);
				getArc(arcId).setCoverage(coverage);
				if(maxConverage<coverage)
					maxConverage=coverage;
				arcId++;
			}
		}

		/*********Set Solution Start Nodes*********/
		readLine(graphFile,line);
		split(splittedLine, line, is_any_of(MPCUtil::SPACE_SEPARATOR));
		for(unsigned j=0;j<splittedLine.size();j++){
			addStartSolNodeId(MPCUtil::getIntValue(splittedLine[j]));
		}

		/*********Set Solution End Nodes*********/
		readLine(graphFile,line);
		split(splittedLine, line, is_any_of(MPCUtil::SPACE_SEPARATOR));
		for(unsigned j=0;j<splittedLine.size();j++){
			addEndSolNodeId(getCopyNodeId(MPCUtil::getIntValue(splittedLine[j])));
		}

		/*********Set subPath Paths*********/
		readLine(graphFile,line);
		int numOfSolPath=MPCUtil::getIntValue(line);
		vector<string> subPathVector;
		for(int i=0;i<numOfSolPath;i++){
			readLine(graphFile,line);
			subPathVector.push_back(line);
		}
		/********Set subpath coverages*********/
		vector<double> subPathCoverages;
		for(int i=0;i<numOfSolPath;i++){
			readLine(graphFile,line);
			subPathCoverages.push_back(MPCUtil::getDoubleValue(line));
		}

		setSubpathConstaints(subPathVector, subPathCoverages);
		graphFile.close();
		MPCLogger::log("Creating graph is finished.");
	}

	// Overloading for the overlaps
	void createMPCGraph(vector<string> &constraints, vector<double> &constraintCoverages, vector<vector<int> > &arcs, vector<int> &sources, vector<int> &sinks) {
		int numOfNodes=constraints.size();
		for (int i=0;i<numOfNodes;i++){
			addNode(i); // create Nodes
			addCopyNode(i);
		}

		/*********Create Arcs*********/
		int arcId=0;
		for(int i=0;i<numOfNodes;i++){
			vector<int> arcVector = arcs.at(i);
			for(unsigned j=0;j<arcVector.size();j++){
				addArc(arcId, getCopyNodeId(i), arcVector.at(j)); // Create arc
				arcId++;
			}
		}
		for(int i=0;i<numOfNodes;i++){
			addArc(arcId, i, getCopyNodeId(i)); // Create arc
			arcId++;
		}

		/*********Set Nodes Costs and lower bounds*********/
		for(int j=0;j<numOfNodes;j++){
			getArc(j,getCopyNodeId(j)).setCost(1.0);
			getArc(j,getCopyNodeId(j)).increaseLowerBound();
		}
		/*********Create the extra node that forces minimum paths *********/
		addNode(numOfNodes);
		addCopyNode(numOfNodes);
		addArc(arcId, numOfNodes, getCopyNodeId(numOfNodes));
		arcId++;

		double sum_of_costs = 0;
		/*********Set Arcs Costs*********/
		arcId=0;
		for(int i=0;i<numOfNodes;i++){
			vector<int> arcVector = arcs.at(i);
			for(unsigned j=0;j<arcVector.size();j++){
				double cost = 0;
				double first_cov = constraintCoverages.at(i);
				double second_cov = constraintCoverages.at(arcVector.at(j));
				cost = abs(1000*(first_cov - second_cov));
				getArc(arcId).setCost(cost);
				arcId++;
				sum_of_costs += cost;
			}
		}
		// Scroll to right place on arcs
		for(int i=0;i<numOfNodes;i++){
			arcId++;
		}

		/********Make the cost of the arc higher than sum of any path to force minimum number of paths *********/
		getArc(arcId).setCost(sum_of_costs);
		arcId++;


		/******* Add arcs from this node to source nodes *******/
		for(unsigned j=0;j<sources.size();j++) {
			addArc(arcId, getCopyNodeId(numOfNodes), sources.at(j));
			arcId++;
		}

		/***The new node is the only source *****/
		addStartSolNodeId(numOfNodes);


		/*********Set Solution End Nodes*********/
		// Constraints containing an end node are end nodes
		for(unsigned j=0;j<sinks.size();j++){
			addEndSolNodeId(getCopyNodeId(sinks.at(j)));
		}

		MPCLogger::log("Creating graph is finished.");
	}

	void setNodesSize(string fileName){
		ifstream nodesFile(fileName.c_str());
		if (!nodesFile.is_open()) { // check for successful opening
			MPCLogger::log("File "+fileName+" does not exist!");
			return ;
		}
		string line;
		vector<string> splittedLine;
		int nodeId=0;
		while(getline(nodesFile,line)){
			split(splittedLine, line, is_any_of(MPCUtil::TAB_SEPARATOR));
			int size=MPCUtil::getIntValue(splittedLine[2])-MPCUtil::getIntValue(splittedLine[1]);
			nodeMap.find(nodeId)->second.setSize(size);
			nodeId++;
		}
	}
	void convertSubpathConstraintIntoEdge(){
		int sourceId,targetId;
		int arcId;
		set<int> subpathNodeIdSet;
		for(vector< vector<int> >::iterator subpathItr=subpathConstraintVector.begin();subpathItr!=subpathConstraintVector.end();subpathItr++){
			vector<int> subpath=*subpathItr;
			for(vector<int>::iterator itr=subpath.begin();itr!=subpath.end();itr++){
				subpathNodeIdSet.insert(*itr);
			}
			if(subpath.size()==1)
				continue;
			sourceId=getCopyNodeId(subpath[0]);
			targetId=subpath[subpath.size()-1];
			arcId=arcMap.size();
			addArc(arcId,sourceId,targetId);
			MPCArc& arc=getArc(arcId);
			arc.increaseLowerBound();
			arc.setSubpathId(distance(subpathConstraintVector.begin(),subpathItr));
		}
		for(map<int,MPCNode>::iterator nodeItr=nodeMap.begin();nodeItr!=nodeMap.end();nodeItr++){
			MPCNode& node=nodeItr->second;
			if(subpathNodeIdSet.find(node.getId())==subpathNodeIdSet.end()){
				MPCArc& dummyArc=getArc(node.getId(), getCopyNodeId(node.getId()));
				dummyArc.increaseLowerBound();
			}
		}
	}
	void setEdgesCost(){
		double cost;
		for(map<int, MPCArc>::iterator arcItr=arcMap.begin();arcItr!=arcMap.end();arcItr++){
			cost=0;
			MPCArc& arc=arcItr->second;
			MPCNode& sourceNode=nodeMap.find(arc.getSource())->second;
			MPCNode& targetNode=nodeMap.find(arc.getTarget())->second;
			cost+=MPCUtil::calcLogCost(sourceNode.getCoverage(), targetNode.getCoverage());
			cost+=MPCUtil::calcLogCost(sourceNode.getCoverage(), arc.getCoverage());
			cost+=MPCUtil::calcLogCost(targetNode.getCoverage(), arc.getCoverage());
//			cost+=MPCUtil::calcLogCost(sourceNode.getCoverage()/maxConverage, targetNode.getCoverage()/maxConverage);
//			cost+=MPCUtil::calcLogCost(sourceNode.getCoverage()/maxConverage, arc.getCoverage()/maxConverage);
//			cost+=MPCUtil::calcLogCost(targetNode.getCoverage()/maxConverage, arc.getCoverage()/maxConverage);
			//still need to add cost of arcs coverage
			arc.setCost(cost);
		}
	}


public:
	static string readLine(ifstream& graphFile, string& line){
		getline(graphFile,line);
		if(ends_with(line, " ")){
			line=line.substr(0,line.size()-1);
		}
		return line;
	}
	static MPCGraph createMPCGraph(string gFileName, string nodesFileName){
		MPCGraph g;
		g.createMPCGraph(gFileName);
//		g.setNodesSize(nodesFileName);
		g.convertSubpathConstraintIntoEdge();
		g.setEdgesCost();
		return g;
	}

	static MPCGraph createMPCGraphFromConstraints(vector<string>& constraints, vector<double>& constraintCoverages, vector<vector<int> > &arcs, vector<int> &sources, vector<int> &sinks) {
		MPCGraph g;
		g.createMPCGraph(constraints, constraintCoverages, arcs, sources, sinks);
		return g;

	}

	void addNode(int nodeId){
		this->nodeMap.insert(pair<int,MPCNode>(nodeId,MPCNode(nodeId)));
	}
	int getCopyNodeId(int nodeId){
		return -nodeId-10;
	}
	void addCopyNode(int nodeId){
		nodeId=getCopyNodeId(nodeId);
		this->nodeMap.insert(pair<int,MPCNode>(nodeId,MPCNode(nodeId)));
	}
	void addArc(int arcId, int source, int destination){
		this->arcMap.insert(pair<int, MPCArc>(arcId,MPCArc(arcId, source, destination)));
		nodeMap.find(source)->second.addOutgoingArc(arcId);
		nodeMap.find(destination)->second.addIncomingArc(arcId);
	}
	MPCNode& getNode(int nodeId){
		return this->nodeMap.find(nodeId)->second;
	}
	MPCArc& getArc(int arcId){
		return this->arcMap.find(arcId)->second;
	}
	MPCArc& getArc(int source, int target){
		vector<int> outArcs=getNode(source).getOutGoingArcs();
		vector<int> inArcs=getNode(target).getIncomingArcs();
		int arcId=-1;
		for(vector<int>::const_iterator outArcsItr=outArcs.begin();outArcsItr!=outArcs.end(); outArcsItr++){
			for(vector<int>::const_iterator inArcsItr=inArcs.begin();inArcsItr!=inArcs.end(); inArcsItr++){
				if(*outArcsItr==*inArcsItr){
					arcId=*outArcsItr;
					break;
				}
			}
			if(arcId!=-1)
				break;
		}
		return getArc(arcId);
	}
	vector<MPCArc> getArc(const vector<int>& arcIdVector) const{
		vector<MPCArc> arcVector;
		for(vector<int>::const_iterator itr=arcIdVector.begin();itr!=arcIdVector.end();itr++){
			arcVector.push_back(this->arcMap.find(*itr)->second);
		}
		return arcVector;
	}
	MPCNode& getSource(MPCArc& arc){
		return this->getNode(arc.getSource());
	}
	MPCNode& getDestination(MPCArc& arc){
		return this->getNode(arc.getTarget());
	}
	const vector<MPCArc> getInArcs(MPCNode& node) const{
		return this->getArc(node.getIncomingArcs());
	}
	const vector<MPCArc> getOutArcs(MPCNode& node) const{
		return this->getArc(node.getOutGoingArcs());
	}
	void addStartSolNodeId(int nodeId){
		this->startSolNodeSet.insert(nodeId);
	}
	void addEndSolNodeId(int nodeId){
		this->endSolNodeSet.insert(nodeId);
	}

	void addSubPathConstraint(vector<int> path){
		this->subpathConstraintVector.push_back(path);
	}

	int countNodes(){
		return this->nodeMap.size();
	}
	int countArcs(){
		return this->arcMap.size();
	}

	const set<int>& getEndSolNodeSet() const {
		return endSolNodeSet;
	}

	const set<int>& getStartSolNodeSet() const {
		return startSolNodeSet;
	}

	const vector< vector<int> >& getSubpathVector() const {
		return subpathConstraintVector;
	}
	string solve(){
		ListDigraph flowNetwork;
		ListDigraph::ArcMap<int64_t> lowerMap(flowNetwork), upperMap(flowNetwork);
		ListDigraph::ArcMap<int64_t> costMap(flowNetwork);
		ListDigraph::NodeMap<int64_t> supplyMap(flowNetwork);

		// adding star source and sink
		ListDigraph::Node s_star = flowNetwork.addNode();
		MPCNode mpcSource(-1);
		mpcSource.setLemonId(flowNetwork.id(s_star));
		flowNodeToMPCNode.insert(pair<int,MPCNode>(flowNetwork.id(s_star),mpcSource));

		ListDigraph::Node t_star = flowNetwork.addNode();
		MPCNode mpcSink(-2);
		mpcSink.setLemonId(flowNetwork.id(t_star));
		flowNodeToMPCNode.insert(pair<int,MPCNode>(flowNetwork.id(t_star),mpcSink));

		ListDigraph::Arc a = flowNetwork.addArc(s_star,t_star); /**/ 
		lowerMap.set(a,0); 
		upperMap.set(a,int64_t_MAX); 
		costMap.set(a,0);
		MPCArc mpcSourceSinkArc(-1,mpcSource.getId(),mpcSink.getId());
		mpcSourceSinkArc.setLemonId(flowNetwork.id(a));
		flowArcToMPCArc.insert(pair<int,MPCArc> (flowNetwork.id(a),mpcSourceSinkArc));
		// Adding the nodes
		for(map<int,MPCNode>::iterator nodeItr=nodeMap.begin();nodeItr!=nodeMap.end();nodeItr++){
			MPCNode& mpcNode=nodeItr->second;
			ListDigraph::Node lemonNode = flowNetwork.addNode(); /**/ 
			supplyMap.set(lemonNode,0);
			int lemonNodeId=flowNetwork.id(lemonNode);
			mpcNode.setLemonId(lemonNodeId); // Setting the node lemon Id
			flowNodeToMPCNode.insert(pair<int,MPCNode>(lemonNodeId,mpcNode));
			if(isStartNode(mpcNode) ){
				MPCLogger::log("Start Node: "+MPCUtil::getStringValue(mpcNode.getId()));
				a = flowNetwork.addArc(s_star,lemonNode); /**/ 
				lowerMap.set(a,0); 
				upperMap.set(a,int64_t_MAX); 
				costMap.set(a,0);
				flowArcToMPCArc.insert(pair<int,MPCArc> (flowNetwork.id(a),MPCArc(-1,mpcSource.getId(),mpcNode.getId())));
			}
			if(isEndNode(mpcNode)){
				MPCLogger::log("End Node: "+MPCUtil::getStringValue(mpcNode.getId()));
				a = flowNetwork.addArc(lemonNode,t_star); /**/ 
				lowerMap.set(a,0); 
				upperMap.set(a,int64_t_MAX); costMap.set(a,0);
				flowArcToMPCArc.insert(pair<int,MPCArc> (flowNetwork.id(a),MPCArc(-2,mpcNode.getId(),mpcSink.getId())));
			}

		}

		// Adding the arcs
		for(map<int,MPCArc>::iterator arcItr=arcMap.begin();arcItr!=arcMap.end();arcItr++){
			MPCArc& mpcArc=arcItr->second;
			ListDigraph::Node source=flowNetwork.nodeFromId(getNode(mpcArc.getSource()).getLemonId());
			ListDigraph::Node destination=flowNetwork.nodeFromId(getNode(mpcArc.getTarget()).getLemonId());
			a=flowNetwork.addArc(source,destination); /**/ 
			lowerMap.set(a,mpcArc.getLowerBound()); 
			upperMap.set(a,mpcArc.getUpperBound()); 
			costMap.set(a,mpcArc.getCost());
			mpcArc.setLemonId(flowNetwork.id(a));
			flowArcToMPCArc.insert(pair<int,MPCArc> (flowNetwork.id(a),mpcArc));
		}
		int totalFlow=100000;
		supplyMap.set(s_star,totalFlow);
		supplyMap.set(t_star,-totalFlow);
		MPCLogger::log("*********Flownetwork Info**********");
		for (ListDigraph::ArcIt arc(flowNetwork); arc != INVALID; ++arc)
		{
			int lowerBound=lowerMap[arc];
			long upperBound=upperMap[arc];
			double cost=costMap[arc];
			int sourceId=flowNodeToMPCNode.find(flowNetwork.id(flowNetwork.source(arc)))->second.getId();
			int targetId=flowNodeToMPCNode.find(flowNetwork.id(flowNetwork.target(arc)))->second.getId();
			MPCLogger::log(MPCUtil::getStringValue(sourceId)+" --> "+MPCUtil::getStringValue(targetId)+"   lowerBound="+MPCUtil::getStringValue(lowerBound)+"      upperBound="+MPCUtil::getStringValue(upperBound) +"  cost="+MPCUtil::getStringValue(cost));
		}
		NetworkSimplex<ListDigraph, int64_t> ns(flowNetwork);
		ns.lowerMap(lowerMap).upperMap(upperMap).costMap(costMap).supplyMap(supplyMap);
		MPCLogger::log("Running the flow engine...");
		NetworkSimplex<ListDigraph, int64_t>::ProblemType res = ns.run(NetworkSimplex<ListDigraph, int64_t>:: CANDIDATE_LIST); //ALTERING_LIST  FIRST_ELIGIBLE BEST_ELIGIBLE BLOCK_SEARCH CANDIDATE_LIST
		if (res != NetworkSimplex<ListDigraph, int64_t>::OPTIMAL)
		{
			MPCLogger::log("ERROR: flow not found");
			return "";
		}
		MPCLogger::log("Finished running the flow engine. ");
		ListDigraph::ArcMap<int64_t> flowMap(flowNetwork); ns.flowMap(flowMap);

		for (ListDigraph::ArcIt arc(flowNetwork); arc != INVALID; ++arc)
		{
			int lowerBound=lowerMap[arc];
			int flow=flowMap[arc];
			double cost=costMap[arc];
			int sourceId=flowNodeToMPCNode.find(flowNetwork.id(flowNetwork.source(arc)))->second.getId();
			int targetId=flowNodeToMPCNode.find(flowNetwork.id(flowNetwork.target(arc)))->second.getId();
			MPCLogger::log(MPCUtil::getStringValue(sourceId)+" --> "+MPCUtil::getStringValue(targetId)
				+"   lowerBound="+MPCUtil::getStringValue(lowerBound)+"   flow="+MPCUtil::getStringValue(flow)+
				"  cost="+MPCUtil::getStringValue(cost));
		}
		string solution="";
		while(totalFlow >0){
			createSolPath(flowNetwork, flowMap, s_star, solution, totalFlow);
		}
		this->solution=solution;
		return solution;
	}

	bool isStartNode(MPCNode& node) const {
		if(node.getIncomingArcs().size()==0 || startSolNodeSet.find(node.getId())!=startSolNodeSet.end())
			return true;
		return false;
	}
	bool isEndNode(MPCNode& node) const {
		if(node.getOutGoingArcs().size()==0 || endSolNodeSet.find(node.getId())!=endSolNodeSet.end())
			return true;
		return false;
	}

	const map<int, MPCArc>& getArcMap() const {
		return arcMap;
	}

	const map<int, MPCNode>& getNodeMap() const {
		return nodeMap;
	}
};

#endif /* MPCGRAPH_H_ */
