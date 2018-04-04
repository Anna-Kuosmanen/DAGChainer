/*
 * SamReader.cpp
 *
 * Created on March 7th 2018
 *	Author: aekuosma
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

//TODO Could add a check that the alignments are clustered by chromosome
// (if meet a chrom that's not lastchrom but is already in the map, then they're not)
bool SamReader::Open(std::string filename) {
	this->samstream.open(filename);

	if(!this->isOpen()) {
		std::cerr << "SamReader error: could not open the specified SAM file." << std::endl;
		return false;
	}

	std::string line;

	std::string lastchrom = "";
	// Saves the start of the line about to be read next
	std::streampos lastpos;

	while(getline(this->samstream, line)) {

		std::vector<std::string> parts = split(line, '\t');
		if(parts.at(2) == lastchrom || parts.at(0)[0] == '@') {
			lastpos = this->samstream.tellg();
			continue;
		}

		chromStarts.insert(std::make_pair(parts.at(2),lastpos));
		lastchrom = parts.at(2);
		lastpos = this->samstream.tellg();

	}
	// Clears EOF
	this->samstream.clear();
	this->samstream.seekg(0, std::ios_base::beg);

	return true;
}

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

