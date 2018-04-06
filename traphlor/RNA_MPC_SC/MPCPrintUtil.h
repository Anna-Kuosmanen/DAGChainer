/*
 * MPCPrintUtil.h
 *
 *  Created on: Feb 14, 2014
 *      Author: ahmedsobih
 */

#ifndef MPCPRINTUTIL_H_
#define MPCPRINTUTIL_H_

class MPCPrintUtil {
public:
	MPCPrintUtil();
	static void printGraph(MPCGraph& g){
		//Print nodes
		MPCLogger::log(MPCUtil::getStringValue(g.countNodes()));
		string line="";
		map<int, MPCNode> nodeMap=g.getNodeMap();
		for(map<int, MPCNode>::const_iterator itr=nodeMap.begin();itr!=nodeMap.end();itr++){
			MPCNode node=itr->second;
			line=MPCUtil::getStringValue(node.getId())+" --> ";
			for(unsigned j=0;j<node.getOutGoingArcs().size();j++){
				line=line+MPCUtil::getStringValue(g.getArc(node.getOutGoingArcs()[j]).getTarget())+" ";
			}
			MPCLogger::log(line);
		}
		//Print nodes end

		//Print nodes coverage
		line="";
		for(map<int, MPCNode>::const_iterator itr=nodeMap.begin();itr!=nodeMap.end();itr++){
			MPCNode node=itr->second;
			if(node.getId()>-1)
				line=line+MPCUtil::getStringValue(node.getCoverage())+" ";
		}
		MPCLogger::log(line);
		//Print nodes coverage end

		//Print arcs coverage

		for(map<int, MPCNode>::const_iterator itr=nodeMap.begin();itr!=nodeMap.end();itr++){
			line="";
			MPCNode node=itr->second;
			for(unsigned j=0;j<node.getOutGoingArcs().size();j++){
				line=line+MPCUtil::getStringValue(g.getArc(node.getOutGoingArcs()[j]).getCoverage())+" ";
			}
			MPCLogger::log(line);
		}
		//Print arcs coverage end

		//Print start solution nodes
		line="";
		for(set<int>::iterator itr=g.getStartSolNodeSet().begin();itr!=g.getStartSolNodeSet().end();itr++){
			line=line+MPCUtil::getStringValue(*itr)+" ";
		}
		MPCLogger::log(line);

		//Print start solution nodes end

		//Print end solution nodes
		line="";
		for(set<int>::iterator itr=g.getEndSolNodeSet().begin();itr!=g.getEndSolNodeSet().end();itr++){
			line=line+MPCUtil::getStringValue(*itr)+" ";
		}
		MPCLogger::log(line);
		//Print end solution nodes end

		//Print subpath
		MPCLogger::log(MPCUtil::getStringValue((int)g.getSubpathVector().size()));
		for(unsigned i=0;i<g.getSubpathVector().size();i++){
			line="";
			vector<int> subPath=g.getSubpathVector()[i];
			for(unsigned j=0;j<subPath.size();j++){
				line=line+MPCUtil::getStringValue(subPath[j])+" ";
			}
			MPCLogger::log(line);
		}

		//Print Print subpath end

		for(map<int, MPCArc>::const_iterator  arcItr=g.getArcMap().begin();arcItr!=g.getArcMap().end();arcItr++){
			MPCArc arc=arcItr->second;
			double cost=arc.getCost();
			int lowerBound=arc.getLowerBound();
			int sourceId=arc.getSource();
			int targetId=arc.getTarget();

			MPCLogger::log(MPCUtil::getStringValue(sourceId)+" --> "+MPCUtil::getStringValue(targetId)+"   lowerBound="+MPCUtil::getStringValue(lowerBound)+"  cost="+MPCUtil::getStringValue(cost));
		}
	}
	static void writeSolutionPath(const string solution, string fileName){
		ofstream solFile((fileName+".sol").c_str());
		MPCLogger::log("The solution paths are: ");
		MPCLogger::log(solution);
		solFile<<solution;
		solFile.close();
	}
};


#endif /* MPCPRINTUTIL_H_ */
