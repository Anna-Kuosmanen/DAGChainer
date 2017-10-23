/*
 * Tuple.h 
 *
 * Created on: Oct 9th 2017
 *		Author: aekuosma
 *
 */

#ifndef TUPLE_H_
#define TUPLE_H_
 
#include <vector>
#include <iostream>
#include <sstream>
#include <cstdlib> 

// Tuples are of form (P,[c,...,d])
class Tuple {

public:	
	// Path
	std::vector<int> P;
	
	int PFirst;
	int PLast;
	
	int c;
	int d;
	
	// Best ordered coverage for some set of tuples that ends in this tuple
	int C;
	// Pointer to the previous tuple (in the chain of best ordered coverage)
	Tuple* previous;
	
	Tuple(std::vector<int> P, int c, int d) {
	
		this->P = P;
		this->c = c;
		this->d = d;
		
		this->PFirst = P.at(0);
		this->PLast = P.at(P.size()-1);

		this->previous = NULL;
		this->C = 0;
	}
	
	std::string toString() {
		std::stringstream sstm;
		sstm << "((";
		for(int i=0;i<P.size()-1;i++)
			sstm << P.at(i) << ", ";
		sstm << this->PLast << "), " << this->c << ", " << this->d << ")";

		return sstm.str();
	}
	bool equals(const Tuple& other) {

		if(this->c != other.c)
			return false;
		if(this->d != other.d)
			return false;

		if(this->P.size() != other.P.size())
			return false;

		for(int i=0;i<this->P.size();i++) {
			if(this->P.at(i) != other.P.at(i)) {
				return false;
			}
		}

		return true;
	}
};

#endif /* TUPLE_H_ */
