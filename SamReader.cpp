/*
 * SamReader.cpp
 *
 * Created on March 7th 2018
 *	Author: Anna Kuosmanen
 *
 */

#include "SamReader.h"
#include "utils.h"

void SamReader::goToChrom(std::string reference) {

	std::streampos pos = chromStarts[reference];

	// Clear EOF
	this->samstream.clear();
	this->samstream.seekg(pos, std::ios_base::beg);

}

SamReader::SamReader() {}

SamReader::~SamReader() {
	if(this->isOpen())
		this->Close();

	this->chromStarts.clear();

}

// Check that the file is ordered, otherwise random access doesn't work
bool SamReader::Open(std::string filename) {
	this->samstream.open(filename);

	if(!this->isOpen()) {
		std::cerr << "SamReader error: could not open the specified SAM file." << std::endl;
		return false;
	}

	bool ordered = true;

	std::string line;

	std::string lastchrom = "";

	// For checking that order is maintained
	std::vector<std::string> seenchroms;
	long lastcoord = 0;

	// Saves the start of the line about to be read next
	std::streampos lastpos;

	while(getline(this->samstream, line)) {

		if(line == "")
			break;

		if(!ordered)
			break;

		std::vector<std::string> parts = split(line, '\t');

		if(parts.at(0)[0] == '@') {
			lastpos = this->samstream.tellg();
			continue;
		}

		if(parts.at(2) == lastchrom) {
			// Check that we're in order within the chromosome
			if(std::atol(parts.at(3).c_str()) >= lastcoord) {
				lastpos = this->samstream.tellg();
				lastcoord = std::atol(parts.at(3).c_str());
				continue;
			}
			else {
				ordered = false;
				break;
			}
		}
		else {
			// New chromosome started, check that we haven't seen it already
			for(unsigned i=0;i<seenchroms.size();i++) {
				if(seenchroms.at(i) == parts.at(2)) {
					ordered = false;
					break;
				}
			}

			chromStarts.insert(std::make_pair(parts.at(2),lastpos));
			lastchrom = parts.at(2);
			seenchroms.push_back(lastchrom);
			lastcoord = std::atol(parts.at(3).c_str());
			lastpos = this->samstream.tellg();
		}

	}
	// Clears EOF
	this->samstream.clear();
	this->samstream.seekg(0, std::ios_base::beg);

	if(ordered)
		return true;
	else {
		this->samstream.close();
		std::cerr << "SamReader error: The SAM file is not sorted." << std::endl;
		return false;
	}
}

// TODO Do some checks that the index values match the file
bool SamReader::Open(std::string filename, std::string index) {
	this->samstream.open(filename);

	std::ifstream indexstream(index);

	if(!this->isOpen()) {
		std::cerr << "SamReader error: could not open the specified SAM file." << std::endl;
		return false;
	}

	if(!indexstream.is_open()) {
		std::cerr << "SamReader error: could not open the specified index file." << std::endl;
		return false;
	}

	std::string line;

	while(getline(indexstream, line)) {

		std::vector<std::string> parts = split(line, '\t');
		chromStarts.insert(std::make_pair(parts.at(0), std::atol(parts.at(1).c_str())));

	}

	indexstream.close();

	return true;
}

// Writes in form "chromName \t position"
bool SamReader::writeIndex(std::string filename) {

	std::ofstream output(filename);

	if(!output.is_open()) {
		std::cerr << "SamReader error: could not open a file for writing the index." << std::endl;
		return false;
	}

	for(auto it = chromStarts.begin(); it != chromStarts.end(); ++it)
		output << it->first << "\t" << it->second << std::endl;

	output.close();

	return true;
}

void SamReader::Close() {
	this->samstream.close();
}

bool SamReader::isOpen() {
	return this->samstream.is_open();
}

SamAlignment SamReader::readAlignment() {

	std::string line;

	if(!this->samstream.eof()) {
		getline(this->samstream, line);
		while(line[0] == '@' && !this->samstream.eof())
			getline(this->samstream, line);
	}
	if(line != "") {
		return SamAlignment(line);
	}
	else {
		return SamAlignment();
	}
}

std::vector<SamAlignment> SamReader::readAlignmentsRegion(std::string reference) {

	this->goToChrom(reference);

	std::vector<SamAlignment> alignments;

	std::streampos lastrowstart = this->samstream.tellg();

	while(true) {
		SamAlignment alignment = this->readAlignment();

		// alignment.isValid() returns false if there was some problem parsing the next alignment
		if(!alignment.isValid())
			break;

		if(alignment.getReferenceName() != reference)
			break;

		alignments.push_back(alignment);
		lastrowstart = this->samstream.tellg();
	}
	// Rewind to the start of first non-eligible alignment
	this->samstream.seekg(lastrowstart, std::ios_base::beg);

	return alignments;

}

std::vector<SamAlignment> SamReader::readAlignmentsRegion(std::string reference, long start, long end) {

	this->goToChrom(reference);

	std::vector<SamAlignment> alignments;

	std::streampos lastrowstart = this->samstream.tellg();

	while(true) {
		SamAlignment alignment = this->readAlignment();

		if(!alignment.isValid())
			break;

		if(alignment.getReferenceName() != reference || alignment.getReferenceStart() > end)
			break;
		if(alignment.getReferenceEnd() >= start)
			alignments.push_back(alignment);

		lastrowstart = this->samstream.tellg();
	}
	// Rewind to the start of first non-eligible alignment
	this->samstream.seekg(lastrowstart, std::ios_base::beg);

	return alignments;


}

std::vector<std::string> SamReader::getReferences() {

	std::unordered_map<std::string, std::streampos>::iterator it;

	std::vector<std::string> chroms;

	for(it=this->chromStarts.begin();it!=this->chromStarts.end();it++)
		chroms.push_back(it->first);

	return chroms;

}

