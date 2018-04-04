/*
 *
 * utils.h
 *
 * Created on: February 19th 2018
 *      Author: aekuosma
 *
 * Implements some utility functions.
 */

#ifndef UTILS_H_
#define UTILS_H_

#include <iostream>
#include <sstream>
#include <vector>

template<typename Out>
static void split(const std::string &s, char delim, Out result) {
        std::stringstream ss;
        ss.str(s);
        std::string item;
        while (std::getline(ss, item, delim)) {
        *(result++) = item;
        }
}

static std::vector<std::string> split(const std::string &s, char delim) {
        std::vector<std::string> elems;
        split(s, delim, std::back_inserter(elems));
        return elems;
}

static std::string reverseComplement(std::string seq) {

	std::string newseq;

	for(int i=seq.size()-1;i>=0;i--) {
		if(seq[i] == 'A')
			newseq += "T";
		else if(seq[i] == 'C')
			newseq += "G";
		else if(seq[i] == 'G')
			newseq += "C";
		else if(seq[i] == 'T')
			newseq += "A";
		else
			newseq += seq[i];

	}

	return newseq;
}

static char complement(char base) {

	if(base=='A' || base=='a')
		return 'T';
	if(base=='C' ||base=='c')
		return 'G';
	if(base=='G' || base=='g')
		return 'C';
	if(base=='T' || base=='t')
		return 'A';
	else
		return base;

}

#endif /* UTILS_H_ */
