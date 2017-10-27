/*
 * Vertex.cpp
 *
 * Created on: Sep 16th 2017
 *	Author: aekuosma
 *
 */
 
#include "Vertex.h"

#include <vector>
#include <cstddef>
#include <algorithm>

Vertex::Vertex(int id, char label, int patlen){
	this->id=id;
	this->label = label;
	for(int i=0;i<patlen;i++){
		this->D.push_back(0);
		this->backpointers.push_back(NULL);
	}
}

int Vertex::getId() {
	return this->id;
}

char Vertex::getLabel() {
	return this->label;
}

// Adds edge from this vertex to v
void Vertex::addEdgeTo(Vertex* v) {
	this->outNeighbors.push_back(v);
	v->inNeighbors.push_back(this);
}

// Adds edge from v to this vertex
void Vertex::addEdgeFrom(Vertex* v) {
	this->inNeighbors.push_back(v);
	v->outNeighbors.push_back(this);
}

void Vertex::updateDValue(int i, int value) {
	if(value >= 0 && i >= 0)
		this->D[i] = value;
}

int Vertex::getDValue(int i) {
	if(i<D.size())
		return this->D[i];
	else
		return -1;
}

std::vector<Vertex*> Vertex::getOutNeighbors(){
	return this->outNeighbors;
}

std::vector<Vertex*> Vertex::getInNeighbors() {
	return this->inNeighbors;
}

std::vector<Vertex*> Vertex::getBackpointers() {
	return this->backpointers;
}

Vertex* Vertex::getBackpointer(int i) {
	if(i<backpointers.size())
		return this->backpointers.at(i);
	else
		return NULL;
}

// Add backpointer for i'th prefix
void Vertex::addBackpointer(int i, Vertex* v) {
	if(i<backpointers.size())
		this->backpointers.at(i) = v;
}

