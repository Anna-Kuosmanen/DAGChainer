/*
 * SamAlignment.h
 *
 * Created on: March 5th 2018
 *	Author: aekuosma
 *
 */

#ifndef SAMALIGNMENT_H_
#define SAMALIGNMENT_H_

#include <vector>
#include <tuple>
#include <map>

class SamAlignment {

private:

	std::string qname;
	int flag;
	std::string rname;
	// First aligned position in the reference
	long reference_start;
	// Last aligned position in the reference
	long reference_end;
	int mapq;
	std::string cigarstring;
	std::vector<std::pair<int,char> > cigartuples;
	// The reference of the mate pair (* for not available, = for same as rname)
	std::string rnext;
	// The position of the mate pair (0 if unavailable)
	long pnext;
	long tlen;
	std::string seq;
	// The quality scores for each base of the read
	std::string qual;
	// Tags such as XS:A for strand of splice alignment
	std::map<std::string, std::string> tags;


public:

	// Dummy, creates an alignment with query name "NULL"
	SamAlignment();

	// Takes a line of SAM file and parses it into alignment
	SamAlignment(std::string description);
	SamAlignment(std::string qname, int flag, std::string rname, long reference_start, long mapq, std::string cigar, std::string rnext, long pnext, long tlen, std::string seq, std::string qual, std::map<std::string, std::string> tags);

	~SamAlignment();

	bool isValid();

	std::string getQueryName();
	void setQueryName(std::string name);

	int getFlag();
	void setFlag(int flag);

	std::string getReferenceName();
	void setReferenceName(std::string name);

	long getReferenceStart();
	void setReferenceStart(long start);

	long getReferenceEnd();
	// Computes reference_end based on reference_start and cigar
	void computeReferenceEnd();
	void setReferenceEnd(long end);

	int getMappingQuality();
	void setMappingQuality(int mapq);

	std::string getCigarString();
	void setCigarString(std::string cigar);

	std::vector<std::pair<int,char> > getCigarTuples();
	void setCigarTuples(std::vector<std::pair<int, char> > newtuples);
	// Splits the cigar string into tuples of (length, operation) and sets the variable cigartuples
	void parseCigarTuples(std::string cigarstring);

	std::string getReferenceForPair();
	void setReferenceForPair(std::string refpair);

	long getPositionForPair();
	void setPositionForPair(long pos);

	long getSignedObservedTemplateLength();
	void setSignedObservedTemplateLength(long len);

	std::string getSequence();
	void setSequence(std::string seq);

	std::string getQuality();
	void setQuality(std::string qual);

	std::map<std::string, std::string> getTags();
	bool hasTag(std::string tag);

	// Do not call this before have checked that the tag is there!
	std::string getValueForTag(std::string tag);
	void setTags(std::map<std::string,std::string> tags);
	void addToTags(std::string tag, std::string value);

	// Various functions to query the flag field
	// In order of the bit, from right to left
	// Is paired
	bool isPaired();
	// Both parts of the pair are aligned
	bool isProperPair();
	bool isUnmapped();
	bool isPairUnmapped();
	bool isReverse();
	bool isPairReverse();
	// Note that isRead1 and isRead2 are only valid for paired, for single reads they're 0
	bool isRead1();
	bool isRead2();
	bool isSecondary();
	bool isQualityControlFail();
	bool isDuplicate();
	bool isSupplementary();

	std::string toString();

};

#endif /* SAMALIGNMENT_H_ */
