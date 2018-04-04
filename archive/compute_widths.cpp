#include "ColinearSolver.h"
#include "BruteForceSolver.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

int main(int argc, char **argv) {

        if(argc < 7) {
                std::cout << "Give the graph file, nodes file, pattern length, anchor length threshold, pattern file and output file prefix." << std::endl;
                return 1;
        }

        std::string graphFile = argv[1];
        std::string nodesFile = argv[2];
        int patlen = std::atoi(argv[3]);

        int threshold = std::atoi(argv[4]);
        std::string thresstring = argv[4];
        std::string patternfile = argv[5];
        std::string outputfile = argv[6];

        std::string brutefile = outputfile + ".brute";
        std::string cofile = outputfile + ".colinear";

        SequenceGraph* SGraph = new SequenceGraph();

        SGraph->readSplicingGraphFile(graphFile, nodesFile, patlen);

        ColinearSolver* cosolver = new ColinearSolver(SGraph);


        std::string tempgraphfile = "flowgraph." + thresstring + "anchor.tmp";

        cosolver->convertSequenceGraphToFlowGraph(tempgraphfile);

        cosolver->solveForPaths(tempgraphfile);


	// Read enough of pattern file to  get the gene ID
        std::ifstream patin(patternfile.c_str());

        std::string line;

        while(getline(patin,line)) {

		break;
	}

	std::string id = line.substr(0,line.size()-2);

	std::cout << id <<  "\t" << cosolver->getPathCoverSize() << std::endl;

}
