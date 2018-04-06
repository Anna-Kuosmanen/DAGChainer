/*
 * MaximalExactMatch.cpp
 *
 * Created on March 20th 2018
 *	Author: aekuosma
 *
 */

#include "MaximalExactMatch.h"



MaximalExactMatch::MaximalExactMatch(int start_, int end_, gcsa::range_type range_) {
	this->range = range_;
	this->begin = start_;
	this->end = end_;

}

