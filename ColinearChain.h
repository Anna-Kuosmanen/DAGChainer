/*
 * ColinearChain.h
 *
 * Created on March 19th 2018
 *	Author: Anna Kuosmanen
 *
 * Basically a vector of tuples and the ordered coverage of the chain.
 *
 */

#ifndef COLINEARCHAIN_H_
#define COLINEARCHAIN_H_

#include "Tuple.h"

#include <vector>

class ColinearChain {

public:

	std::vector<Tuple*> chain;
	int coverageScore;

	ColinearChain(std::vector<Tuple*> chain_, int coverageScore_) {
		this->chain = chain_;
		this->coverageScore = coverageScore_;
	}

	// Dummy to use as "no chain found"
	ColinearChain() {
		this->coverageScore = -1;

	}

	std::string toString() {

		std::stringstream sstm;

		sstm << "Score: " << this->coverageScore << ". Chain: ";

		for(unsigned i=0;i<chain.size();i++)
			sstm << chain.at(i)->toString() << ",";

		return sstm.str();

	}

};

#endif /* COLINEAR_CHAIN_H_ */
