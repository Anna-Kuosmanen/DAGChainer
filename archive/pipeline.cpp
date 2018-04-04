#include "ColinearSolver.h"
#include "BruteForceSolver.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <chrono>

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

	bool readsuccess = SGraph->readSplicingGraphFile(graphFile, nodesFile, patlen);

	if(!readsuccess) {
		std::cerr << "Failed to create the sequence graph, check the existance of the graph file and nodes file." << std::endl;
		return 0;

	}

	ColinearSolver* cosolver = new ColinearSolver(SGraph);
	BruteForceSolver *brutesolver = new BruteForceSolver(SGraph);

	std::string tempgraphfile = "flowgraph." + thresstring + "anchor.tmp";

	cosolver->convertSequenceGraphToFlowGraph(tempgraphfile);

	std::cerr << "Case\tTime\t|V|\tk\tN" << std::endl;


	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	cosolver->solveForPaths(tempgraphfile);

	std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();

	std::cerr << "Solve_path\t" <<  std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "\t" << SGraph->getNoOfVertices() << "\t" << cosolver->getPathCoverSize() << "\t-" << std::endl;
	begin = std::chrono::steady_clock::now();

	cosolver->computeForward();

        end= std::chrono::steady_clock::now();

        std::cerr << "Compute_forward\t" <<  std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "\t" << SGraph->getNoOfVertices() << "\t" << cosolver->getPathCoverSize() << "\t-" << std::endl;

	std::ifstream patin(patternfile.c_str());
	std::ofstream outco(cofile.c_str());

	if(!patin.good()) {
		std::cerr << "Failed to open the pattern file, exiting." << std::endl;
		return 0;

	}

	if(!outco.good()) {
		std::cerr << "Failed to open the output file, exiting." << std::endl;
		return 0;
	}

	std::ofstream outbrute(brutefile.c_str());

	outco << "##Vertices:\t" << SGraph->getNoOfVertices() << std::endl;
	outco << "##Paths:\t" << cosolver->getPathCoverSize() << std::endl;

	outbrute << "##Vertices:\t" << SGraph->getNoOfVertices() << std::endl;

	std::string line;

	while(getline(patin,line)) {

		std::string id = line;

		getline(patin, line);


        	begin = std::chrono::steady_clock::now();

		std::vector<Tuple*> anchors;

		SGraph->findAnchors(anchors, line, threshold);

		end= std::chrono::steady_clock::now();
		std::cerr << "Find_anchors/" << id << "\t" <<  std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "\t" << SGraph->getNoOfVertices() << "\t" << cosolver->getPathCoverSize() << "\t" << anchors.size()  << std::endl;


		outco << id << "#anchors=" << anchors.size() << std::endl;
		outbrute << id << "#anchors=" << anchors.size() << std::endl;


		begin = std::chrono::steady_clock::now();

		std::vector<Tuple*> cochain;

		 cosolver->solveForAnchors(anchors, cochain);

		end= std::chrono::steady_clock::now();
		std::cerr << "Fancy_algorithm/" << id << "\t" <<  std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "\t" << SGraph->getNoOfVertices() << "\t" << cosolver->getPathCoverSize() << "\t" << anchors.size()  << std::endl;

		begin = std::chrono::steady_clock::now();

		std::vector<Tuple*> brutechain;
		brutesolver->solveForAnchors(anchors, brutechain);

		end= std::chrono::steady_clock::now();
		std::cerr << "Naive_algorithm/" << id << "\t" <<  std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "\t" << SGraph->getNoOfVertices() << "\t" << cosolver->getPathCoverSize() << "\t" << anchors.size()  << std::endl;

	

		for(unsigned i=0;i<cochain.size();i++)
			outco << cochain.at(i)->toString() << ", ";
		outco << std::endl;

		for(unsigned i=0;i<brutechain.size();i++)
			outbrute << brutechain.at(i)->toString() << ", ";
		outbrute << std::endl; 


		// Clean-up the Tuples (chains use subset of anchors, this is enough)
		for(unsigned i=0;i<anchors.size();i++)
			delete anchors.at(i);
	}


	delete SGraph;
	delete brutesolver;
	delete cosolver;


	patin.close();
	outco.close();
	outbrute.close();

}

