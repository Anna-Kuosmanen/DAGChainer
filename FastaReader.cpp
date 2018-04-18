/*
 * FastaReader.cpp
 *
 * Created March 18th 2018
 *	Author: Anna Kuosmanen
 *
 */

#include <algorithm>

#include "FastaReader.h"


FastaReader::FastaReader() {}

FastaReader::~FastaReader() {

	if(this->fastastream.is_open())
		this->Close();

}

bool FastaReader::Open(std::string filename) {

	this->fastastream.open(filename.c_str());

	if(!this->fastastream.is_open()) {
		std::cerr << "FastaReader error: failed to open the input file." << std::endl;
		return false;
	}

	return true;
}

void FastaReader::Close() {
	this->fastastream.close();
}

FastaEntry FastaReader::next() {

	std::string line = "";
	std::string id = "";
	std::string seq = "";

	std::streampos lastrowstart = this->fastastream.tellg();

	// The sequence can span multiple lines
	while(!this->fastastream.eof()) {
		getline(this->fastastream,line);
		if(line == "")
			break;
		// ID line
		if(line.at(0) == '>') {
			// If id is empty, this is the id
			if(id == "") {
				id = line.substr(1, std::string::npos);
				lastrowstart = this->fastastream.tellg();
			}
			// Else we read beyond interest, break (rewind is done at end)
			else {
				break;
			}

		}
		else {
			seq += line;
			lastrowstart = this->fastastream.tellg();
		}
	}

	// Uppercase in case there are lowercase bases, those cause issues
	std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);


	this->fastastream.seekg(lastrowstart, std::ios_base::beg);

	if(seq != "") {
		return FastaEntry(id, seq);
	}

	// Return entry with empty id and seq if there was nothing to read
	else {
		return FastaEntry();
	}


}
