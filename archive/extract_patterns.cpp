// Input: graph file, nodes file, number and length of patterns.
// Output: patterns printed to standard output


#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <stdlib.h>


#include "SequenceGraph.h"


int main(int argc, char **argv) {

	if(argc < 6) {
		std::cout << "Give the graph file, nodes file, pattern name prefix, number and length of patterns." << std::endl;
		return 1;
	}

	std::string graphFile = argv[1];
	std::string nodesFile = argv[2];
	// This isn't used, keeping it for the sake of arguments being the same as in pipelines
	std::string patprefix = argv[3];
	int noOfPatterns = std::atoi(argv[4]);
	int length = std::atoi(argv[5]);


	SequenceGraph* SGraph = new SequenceGraph();
	// Choose the initializing patlen as 1 for speed, won't need it
	bool readsuccess = SGraph->readSplicingGraphFile(graphFile, nodesFile, 1);

	if(!readsuccess) {
		std::cerr << "Failed to create the sequence graph, check the existance of the graph file and node file." << std::endl;
		return 0;
	}

	std::vector<std::string> patterns;

	// Call it quits if it's turning out to be too hard to find paths of this length
	int maxattempts = 100*noOfPatterns;

	int attempts = 0;

	while(patterns.size() < noOfPatterns) {

		if(attempts > maxattempts) {
			std::cout << "Tried too many times, this isn't working!" << std::endl;
			return 0;
		}

		attempts++;

		// Starting farther than length of bases from end is pointless
		int choice = rand() % (SGraph->getNoOfVertices()-length);

		Vertex* current = SGraph->getVertex(choice);

		std::string pattern;

		pattern += current->getLabel();

		while(pattern.size() < length && current->getOutNeighbors().size() > 0) {
			std::vector<Vertex*> neighbors = current->getOutNeighbors();

			int neighborchoice = rand() % (neighbors.size());

			current = neighbors.at(neighborchoice);

			pattern += current->getLabel();


		}

		if(pattern.size() >= length)
			patterns.push_back(pattern);

	}


	// Print
	for(int i=0;i<noOfPatterns;i++) {
		std::cout << ">" << patprefix << "." << i << std::endl;
		std::cout << patterns.at(i) << std::endl;
	}

}
