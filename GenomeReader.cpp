/*
 * GenomeReader.cpp
 *
 * Created on March 14th 2018
 *	Author: Anna Kuosmanen
 *
 */

#include <algorithm>

#include "GenomeReader.h"
#include "utils.h"

GenomeReader::GenomeReader() {

}

GenomeReader::~GenomeReader() {

	if(this->isOpen())
		this->Close();
	this->chromStarts.clear();

}

bool GenomeReader::Open(std::string filename) {

	this->genomestream.open(filename);

	if(!this->isOpen()) {
		std::cerr << "GenomeReader error: could not open the input file." << std::endl;
		return false;
	}

	std::string line;

	// To check that all lines comply to length;
	int lastseqlength = -1;
	// Makes sure only the last line of an entry can be different length
	bool lastsegmentline = false;

	std::streampos lastpos = this->genomestream.tellg();

	while(getline(this->genomestream, line)) {

		// ID line, save start to map
		if(line[0] == '>') {
			this->chromStarts.insert(std::make_pair(line.substr(1,std::string::npos), this->genomestream.tellg()));
			this->chroms.push_back(line.substr(1,std::string::npos));
			lastsegmentline = false;
			lastpos = this->genomestream.tellg();
		}
		// Sequence line, check that the lengths match
		else {
			// lastsegmentline flags if the length of the line was different from before, and resets when ID is found
			if(lastsegmentline) {
				std::cerr << "GenomeReader error: rows of the sequence are not of same length." << std::endl;
					std::cerr << "Random access requires all lines to be of same length (except for the last line of each sequence)." << std::endl;
				return false;
			}
			// This is the first sequence line, save it to class variable
			if(lastseqlength == -1) {
				lastseqlength = int(line.size());
				this->seqlinelength = lastseqlength;
				this->pointerlinelength = int(this->genomestream.tellg() - lastpos);
			}
			// If it's different size, mark this as the last segment
			else if(int(line.size()) != this->seqlinelength) {
				lastsegmentline = true;
			}
			else {
				// If the sequence length matches but streampos values don't, there's something up with the newlines
				if(int(this->genomestream.tellg()-lastpos) != this->pointerlinelength) {
					std::cerr << "GenomeReader error: check the newlines in your file (possible mix of Unix and Windows coding)." << std::endl;
					std::cerr << "Random access requires all lines to be of same length." << std::endl;
					return false;
				}
			}
			lastpos = this->genomestream.tellg();
		}
	}
	// Clear EOF
	this->genomestream.clear();

	// Save the end (needed to check that searches in the last chromosome don't go over the end of file)
	this->genomestream.seekg(0, std::ios_base::end);
	this->chromStarts.insert(std::make_pair("END",this->genomestream.tellg()));
	this->chroms.push_back("END");

	this->genomestream.seekg(0, std::ios_base::beg);

	return true;
}

bool GenomeReader::Open(std::string filename, std::string index) {

	this->genomestream.open(filename);

	if(!this->isOpen()) {
		std::cerr << "GenomeReader error: could not open the input file." << std::endl;
		return false;
	}

	std::ifstream indexstream(index);

	if(!indexstream.is_open()) {

		std::cerr << "GenomeReader error: could not open the specified index file." << std::endl;
		return false;

	}

	std::string line;

	getline(indexstream, line);

	if(line[0] != '@') {
		std::cerr << "GenomeReader error: the index is not valid! The first line should contain @ and the length of each line of sequence (excluding newlines)" << std::endl;
		indexstream.close();
		return false;
	}

	seqlinelength = atoi(line.substr(1,std::string::npos).c_str());

	getline(indexstream, line);

	if(line[0] != '@') {
		std::cerr << "GenomeReader error: the index is not valid! The second line should contain @ and the length of each line of sequence (including newlines)" << std::endl;
		indexstream.close();
		return false;
	}

	pointerlinelength = atoi(line.substr(1,std::string::npos).c_str());

	while(getline(indexstream, line)) {

		std::vector<std::string> parts = split(line, '\t');
		this->chromStarts.insert(std::make_pair(parts.at(0), std::atol(parts.at(1).c_str())));
		this->chroms.push_back(parts.at(0));

	}

	indexstream.close();

	return true;
}

// Writes in form "chromName \t position"
bool GenomeReader::writeIndex(std::string filename) {

	std::ofstream output(filename);

	if(!output.is_open()) {
		std::cerr << "GenomeReader error: could not open a file for writing the index." << std::endl;
		return false;
	}

	output << "@" << this->seqlinelength << std::endl;
	output << "@" << this->pointerlinelength << std::endl;

	for(unsigned i=0;i<this->chroms.size();i++) {
		auto it = this->chromStarts.find(chroms.at(i));
		output << it->first << "\t" << it->second << std::endl;
	}

	output.close();
	return true;
}


void GenomeReader::Close() {
	this->genomestream.close();
}

bool GenomeReader::isOpen() {

	return this->genomestream.is_open();

}

bool GenomeReader::hasChromosome(std::string chrom) {

	if(this->chromStarts.find(chrom) != this->chromStarts.end())
		return true;
	else
		return false;
}


std::string GenomeReader::readSequence(std::string chrom, long start, long end) {
	long length = end-start+1;

	if(this->hasChromosome(chrom)) {

		std::streampos startpos = this->chromStarts.find(chrom)->second;
		long diff = (start/this->seqlinelength)*this->pointerlinelength+(start%this->seqlinelength);
		// Check that the position is still within this chromosome
		// TODO This is off by the length of the chromosome name, fix
		auto it = this->chromStarts.find(chrom);
		it++;
		std::streampos nextchromstart = it->second;
		if(startpos+diff+length < nextchromstart) {
			this->genomestream.seekg(startpos+diff, std::ios_base::beg);

			// Read lines till get "length" long sequence
			std::string line;

			std::string seq = "";

			while(long(seq.size()) < length) {
				getline(genomestream,line);
				seq += line;
			}

			seq = seq.substr(0, length);


			return seq;
		}
		else
			return "";
	}
	
	return "";
}

std::string GenomeReader::readSequenceUpperCase(std::string chrom, long start, long end) {

	std::string seq = this->readSequence(chrom, start, end);

	std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);

	return seq;

}


