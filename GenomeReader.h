/*
 * GenomeReader.h
 *
 * Created on March 14th 2018
 *	Author: aekuosma
 *
 * Special case of FASTA reader, indexes genome FASTA files to allow for fast sequence retrieval
 * (without reading the whole genome to memory)
 *
 */

#ifndef GENOMEREADER_H_
#define GENOMEREADER_H_

#include <iostream>
#include <fstream>
#include <unordered_map>

class GenomeReader {

private:

	std::ifstream genomestream;

	std::unordered_map<std::string, std::streampos> chromStarts;
	// Saves the chrom names in order since map doesn't keep order
	std::vector<std::string> chroms;

	// Length of each line of sequence (except the last for each chromosome), excludes newlines
	int seqlinelength;

	// Length of each line of sequence, includes newlines (newline can be either "\n" or "\r\n" depending on system, so this is portable)
	int pointerlinelength;


public:

	GenomeReader();
	~GenomeReader();

	// Slower, creates the index for random access
	bool Open(std::string filename);

	// Faster, reads the index from file
	bool Open(std::string filename, std::string index);

	bool writeIndex(std::string filename);

	void Close();

	bool isOpen();

	// Checks whether the index has given chromosome
	bool hasChromosome(std::string chrom);

	// 0-based, so if using with 1-based coordinates, substract 1
	std::string readSequence(std::string chrom, long start, long end);

	// The same as above, except converts the whole thing into uppercase
	std::string readSequenceUpperCase(std::string chrom, long start, long end);


};

#endif /* GENOMEREADER_H_ */
