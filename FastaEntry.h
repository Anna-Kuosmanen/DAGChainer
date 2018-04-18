/*
 * FastaEntry.h
 * Created on March 18th 2018
 *	Author: Anna Kuosmanen
 *
 * A single FASTA entry for the use of FastaReader.
 *
 */

#ifndef FASTAENTRY_H_
#define FASTAENTRY_H_

#include <string>

class FastaEntry {

public:

        std::string id;
        std::string seq;

	FastaEntry(std::string id_, std::string seq_) {
		this->id = id_;
		this->seq = seq_;
	}

	FastaEntry() {
		this->id = "";
		this->seq = "";
	}
};

#endif /* FASTAENTRY_H_ */
