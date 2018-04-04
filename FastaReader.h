/*
 * FastaReader.h
 *
 * Created on March 18th 2018
 *	Author: aekuosma
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

	FastaEntry next();




};
#endif /* FASTAREADER_H_ */
