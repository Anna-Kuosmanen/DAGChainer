/*
 * Vertex.h
 *
 * Created on: Sep 16th 2017
 *	Author: Anna Kuosmanen
 *
 * This is basically a linked list structure
 * 
 */

#ifndef VERTEX_H_
#define VERTEX_H_

#include <vector>

class Vertex {

private:
	int id;
	char label;
	std::vector<int> D; // The vector holding the lengths of suffix of prefix (for naive anchor finding)
	std::vector<Vertex*> backpointers; // Where the path corresponding to i'th suffix came from (again for naive anchor finding)
	std::vector<Vertex*> outNeighbors;
	std::vector<Vertex*> inNeighbors;


public:
	Vertex(int id, char label, int patlen);

	int getId();

	char getLabel();

	// Adds edge from this vertex to v
	void addEdgeTo(Vertex* v);

	// Adds edge from v to this vertex
	void addEdgeFrom(Vertex* v);

	void updateDValue(int i, int value);

	int getDValue(int i);

	std::vector<Vertex*> getOutNeighbors();

	std::vector<Vertex*> getInNeighbors();

	std::vector<Vertex*> getBackpointers();

	Vertex* getBackpointer(int i);

	// Add backpointer for i'th prefix
	void addBackpointer(int i, Vertex* v);

};
#endif /* VERTEX_H_ */
