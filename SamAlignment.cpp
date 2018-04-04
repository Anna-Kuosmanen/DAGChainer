/*
 * SamAlignment.cpp
 *
 * Created on: March 5th 2018
 *	Author: aekuosma
 *
 */

#include "SamAlignment.h"
#include "utils.h"
#include <iostream>
#include <sstream>
#include <stdlib.h> /* atoi */
#include <ctype.h> /* isdigit */

//Note that positions are saved in 0-based, but printed in 1-based to comply with SAM standard

SamAlignment::SamAlignment() {
	this->setQueryName("NULL");
}

SamAlignment::SamAlignment(std::string description){
	std::vector<std::string> parts = split(description, '\t');

	// Check that the entry at least seems valid (strings for strings, numbers for numbers, etc.)
	bool faulty = false;

	if(parts.size() < 11) {
		faulty = true;
	}
	else {
		// Some field is blank
		for(unsigned i=0;i<=11;i++)
			if(parts.at(i) == "")
				faulty = true;
	}

	// Parts 1, 3, 4, 7 and 8 should be numbers (0-based, query name is 0th part)
	for(unsigned i=0;i<parts.at(1).size();i++)
		if(!isdigit(parts.at(1).at(i)))
			faulty = true;
	for(unsigned i=0;i<parts.at(3).size();i++)
		if(!isdigit(parts.at(3).at(i)))
			faulty = true;
	for(unsigned i=0;i<parts.at(4).size();i++)
		if(!isdigit(parts.at(4).at(i)))
			faulty = true;
	for(unsigned i=0;i<parts.at(7).size();i++)
		if(!isdigit(parts.at(7).at(i)))
			faulty = true;
	for(unsigned i=0;i<parts.at(8).size();i++)
		if(!isdigit(parts.at(8).at(i)))
			faulty = true;


	if(faulty)
		this->setQueryName("NULL");
	else {

		this->setQueryName(parts.at(0));
		this->setFlag(atoi(parts.at(1).c_str()));
		this->setReferenceName(parts.at(2));
		this->setReferenceStart(atol(parts.at(3).c_str())-1);
		this->setMappingQuality(atoi(parts.at(4).c_str()));
		this->setCigarString(parts.at(5));
		this->parseCigarTuples(parts.at(5));
		this->setReferenceForPair(parts.at(6));
		this->setPositionForPair(atol(parts.at(7).c_str())-1);
		this->setSignedObservedTemplateLength(atol(parts.at(8).c_str()));
		this->setSequence(parts.at(9));
		this->setQuality(parts.at(10));

		// Sets reference_end based on reference_start and cigar
		this->computeReferenceEnd();

		// If there are more (not empty space, the line can end with tab),
		// then there are tags
		if(parts.size() > 11 && parts.at(11) != "") {

			std::stringstream sstm;

			for(unsigned i= 11; i<parts.size();i++) {
				if(parts.at(i) == "")
					break;
				// Split on the last ":"
				std::vector<std::string> tagparts = split(parts.at(i), ':');

				for(unsigned t=0;t<tagparts.size()-2;t++)
					sstm << tagparts.at(t) << ":";

				sstm << tagparts.at(tagparts.size()-2);

				std::string key = sstm.str();
				sstm.str("");

				std::string value = tagparts.at(tagparts.size()-1);		

				this->addToTags(key, value);	

			}
		}
	}	
}


// This expects 0-based
SamAlignment::SamAlignment(std::string qname, int flag, std::string rname, long reference_start, long mapq, std::string cigar, std::string rnext, long pnext, long tlen, std::string seq, std::string qual, std::map<std::string,std::string> tags){

	this->setQueryName(qname);
	this->setFlag(flag);
	this->setReferenceName(rname);
	this->setReferenceStart(reference_start);
	this->setMappingQuality(mapq);
	this->setCigarString(cigar);
	this->parseCigarTuples(cigar);
	this->setReferenceForPair(rnext);
	this->setPositionForPair(pnext);
	this->setSignedObservedTemplateLength(tlen);
	this->setSequence(seq);
	this->setQuality(qual);
	this->setTags(tags);
}

SamAlignment::~SamAlignment() {

	this->tags.clear();
	this->cigartuples.clear();

}

bool SamAlignment::isValid() {

	if(this->getQueryName() == "NULL")
		return false;
	else
		return true;

}

std::string SamAlignment::getQueryName() {
	return this->qname;
}

void SamAlignment::setQueryName(std::string name){
	this->qname = name;
}

int SamAlignment::getFlag() {
	return this->flag;
}

void SamAlignment::setFlag(int flag) {
	if(flag >= 0)
		this->flag = flag;
}

std::string SamAlignment::getReferenceName(){
	return this->rname;
}

void SamAlignment::setReferenceName(std::string name){
	this->rname = name;
}

long SamAlignment::getReferenceStart(){
	return this->reference_start;
}

void SamAlignment::setReferenceStart(long start){
	if(start >= 0)
		this->reference_start = start;
}

long SamAlignment::getReferenceEnd(){
	return this->reference_end;
}

// Depending on the convention, reference_end can be smaller than reference_start on the reverse strand, so can't demand that reference_end > reference_start
void SamAlignment::setReferenceEnd(long end){
	if(end >= 0)
		this->reference_end = end;
}

void SamAlignment::computeReferenceEnd() {

	long pos = this->getReferenceStart();

	for(unsigned i=0;i<this->cigartuples.size();i++) {
		std::pair<int,char> part = this->cigartuples.at(i);

		// Inserts, soft and hard clips don't move on the reference, the rest do
		if(part.second == 'I' || part.second == 'S' || part.second == 'H')
			continue;

		pos += part.first;

	}
	// -1 because the last base is at start+length-1
	this->setReferenceEnd(pos-1);

}

int SamAlignment::getMappingQuality(){
	return this->mapq;
}

// Mapping quality is 0-255
void SamAlignment::setMappingQuality(int mapq){
	if(mapq >=0 && mapq <=255)
		this->mapq = mapq;
}

std::string SamAlignment::getCigarString(){
	return this->cigarstring;
}

void SamAlignment::setCigarString(std::string cigar){
	this->cigarstring = cigar;
}

std::vector<std::pair<int,char> > SamAlignment::getCigarTuples(){
	return this->cigartuples;
}

void SamAlignment::setCigarTuples(std::vector<std::pair<int, char> > newtuples){
	this->cigartuples = newtuples;
}

void SamAlignment::parseCigarTuples(std::string cigarstring){

	std::stringstream is(cigarstring);
	char n;

	std::string current = "";

	while(is >> n) {
		if(isdigit(n)) {
			current += n;
		}
		else {
			this->cigartuples.push_back(std::make_pair(atoi(current.c_str()),n));
			current = "";
		}
	}
}

std::string SamAlignment::getReferenceForPair(){
	return this->rnext;
}
void SamAlignment::setReferenceForPair(std::string refpair){
	this->rnext = refpair;
}

long SamAlignment::getPositionForPair(){
	return this->pnext;
}
void SamAlignment::setPositionForPair(long pos){
	if(pos >=0)
		this->pnext = pos;
}

long SamAlignment::getSignedObservedTemplateLength(){
	return this->tlen;
}
void SamAlignment::setSignedObservedTemplateLength(long len){
	this->tlen = len;
}

std::string SamAlignment::getSequence(){
	return this->seq;
}

void SamAlignment::setSequence(std::string seq){
	this->seq = seq;
}

std::string SamAlignment::getQuality(){
	return this->qual;
}

void SamAlignment::setQuality(std::string qual){
	this->qual = qual;
}

std::map<std::string,std::string> SamAlignment::getTags(){
	return this->tags;
}

bool SamAlignment::hasTag(std::string tag) {
	return((this->tags.count(tag))>0);
}

std::string SamAlignment::getValueForTag(std::string tag) {
	return (this->tags.find(tag)->second);
}

void SamAlignment::setTags(std::map<std::string,std::string>  tags){
	this->tags = tags;
}

void SamAlignment::addToTags(std::string tag, std::string value){
	this->tags.insert(std::make_pair(tag, value));
}


bool SamAlignment::isPaired(){
	if((this->flag)%2 == 1)
		return true;
	else
		return false;	
}

bool SamAlignment::isProperPair(){
	if(((this->flag)>>1)%2 == 1)
		return true;
	else
		return false;
}

bool SamAlignment::isUnmapped(){
	if(((this->flag)>>2)%2 == 1)
		return true;
	else
		return false;
}

bool SamAlignment::isPairUnmapped() {
	if(((this->flag)>>3)%2 == 1)
		return true;
	else
		return false;
}

bool SamAlignment::isReverse(){
	if(((this->flag)>>4)%2 == 1)
		return true;
	else
		return false;
}

bool SamAlignment::isPairReverse() {
	if(((this->flag)>>5)%2 == 1)
		return true;
	else
		return false;
}

bool SamAlignment::isRead1(){
	if(((this->flag)>>6)%2 == 1)
		return true;
	else
		return false;
}
bool SamAlignment::isRead2(){
	if(((this->flag)>>7)%2 == 1)
		return true;
	else
		return false;
}

bool SamAlignment::isSecondary(){
	if(((this->flag)>>8)%2 == 1)
		return true;
	else
		return false;

}

bool SamAlignment::isQualityControlFail(){
	if(((this->flag)>>9)%2 == 1)
		return true;
	else
		return false;
}

bool SamAlignment::isDuplicate(){
	if(((this->flag)>>10)%2 == 1)
		return true;
	else
		return false;
}

bool SamAlignment::isSupplementary(){
	if(((this->flag)>>11)%2 == 1)
		return true;
	else
		return false;
}

std::string SamAlignment::toString() {

	std::stringstream sstm;

	sstm << this->getQueryName() << "\t" << this->getFlag() << "\t" << this->getReferenceName() << "\t" << this->getReferenceStart()+1 << "\t" << this-> getMappingQuality() << "\t" << this->getCigarString() << "\t" << this->getReferenceForPair() << "\t" << this->getPositionForPair() << "\t" << this->getSignedObservedTemplateLength() << "\t" << this->getSequence() << "\t" << this->getQuality() << "\t";

	std::map<std::string, std::string> tags = this->getTags();
	std::map<std::string, std::string>::iterator it;

	for(it=tags.begin();it!=tags.end();it++)
		sstm << it->first << ":" << it->second << "\t";

	return sstm.str();
}

