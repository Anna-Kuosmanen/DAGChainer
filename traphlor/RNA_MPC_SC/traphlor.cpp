/*
 * main.cpp
 *
 *  Created on: Feb 11, 2014
 *      Author: ahmedsobih
 */
#ifndef MAIN_CPP_
#define MAIN_CPP_

#include "MPCGraph.h"
#include "MPCUtil.h"
#include "MPCPrintUtil.h"

int main(int argc, char **argv) {
	if (argc<3) {
		MPCLogger::log("Parameters: [MPC_Graph_File] [Nodes_File]" );
		return 0;
	}
	string graphFile=argv[1];
	string nodeFile=argv[2];
	MPCGraph g=MPCGraph::createMPCGraph(graphFile, nodeFile);
	MPCPrintUtil::printGraph(g);
	string solution=g.solve();
	MPCPrintUtil::writeSolutionPath(solution, graphFile);

}
#endif /* MAIN_CPP_ */
