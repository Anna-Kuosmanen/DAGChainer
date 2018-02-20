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




template<typename Out>
void SequenceGraph::split(const std::string &s, char delim, Out result) {
        std::stringstream ss;
        ss.str(s);
        std::string item;
        while (std::getline(ss, item, delim)) {
        *(result++) = item;
        }
}

std::vector<std::string> SequenceGraph::split(const std::string &s, char delim) {
        std::vector<std::string> elems;
        split(s, delim, std::back_inserter(elems));
        return elems;
}


#endif /* UTILS_H_ */
