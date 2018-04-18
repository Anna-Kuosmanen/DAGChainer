/*
 * FastaReader.h
 *
 * Created on March 18th 2018
 *	Author: Anna Kuosmanen
 *
 * Reads FASTA files. Main use is for when the sequence of an entry is split on several rows.
 *
 */

#ifndef FASTAREADER_H_
#define FASTAREADER_H_

#include "FastaEntry.h"

#include <iostream>
#include <fstream>

class FastaReader {


private:

	std::ifstream fastastream;

public:

	FastaReader();
	~FastaReader();

	bool Open(std::string filename);
	void Close();

	// Reads the next FASTA entry. Returns an entry with id and seq of "" if nothing to read anymore
	// The sequence is converted to uppercase if there are any lowercase entries
	FastaEntry next();




};
#endif /* FASTAREADER_H_ */
