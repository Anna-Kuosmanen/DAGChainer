/*
 * MaximalExactMatch.h
 *
 * Created on March 20th 2018
 *	Author: Anna Kuosmanen
 *
 * A struct for a MEM.
 */

#ifndef MAXIMALEXACTMATCH_H_
#define MAXIMALEXACTMATCH_H_

#include "Tuple.h"
#include "gcsa/gcsa.h"
#include "gcsa/lcp.h"

#include <vector>

class MaximalExactMatch {

public:
	// The gcsa range for the MEM
	gcsa::range_type range;
	std::vector<gcsa::node_type> nodes;

	// The start and end of the MEM in the sequence
	int begin;
	int end;

	MaximalExactMatch(int start_, int end_, gcsa::range_type range_);



};

#endif /* MAXIMALEXACTMATCH_H_ */
