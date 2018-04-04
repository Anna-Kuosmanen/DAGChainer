/*
 * SamReader.h
 *
 * Created on March 7th 2018
 *	Author: aekuosma
 *
 */

#ifndef SAMREADER_H_
#define SAMREADER_H_

#include "SamAlignment.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>


class SamReader {

private:

	std::ifstream samstream;

	std::unordered_map<std::string, std::streampos> chromStarts;
	

public:

	SamReader();
	~SamReader();

	// Slower, since reads through the file to create the chromStarts map
	bool Open(std::string filename);

	// Reads the chromStarts from file "index"
	bool Open(std::string filename, std::string index);

	bool writeIndex(std::string filename);

	void Close();

	bool isOpen();

	// Makes region location faster, moves the read pointer to the first read aligned to this chromosome
	void goToChrom(std::string reference);

	SamAlignment readAlignment();

	std::vector<SamAlignment> readAlignmentsRegion(std::string reference);

	std::vector<SamAlignment> readAlignmentsRegion(std::string reference, long start, long end);

	// Lists the chromosomes in the SAM file
	std::vector<std::string> getReferences();



};

#endif /* SAMREADER_H_ */
