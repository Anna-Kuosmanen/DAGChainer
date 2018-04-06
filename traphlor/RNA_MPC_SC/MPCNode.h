/*
 * MPCNode.h
 *
 *  Created on: Feb 11, 2014
 *      Author: ahmedsobih
 */

#ifndef MPCNODE_H_
#define MPCNODE_H_

class MPCNode {
private:
	int id;
	int lemonId;
	double coverage;
	vector<int> incomingArcs;
	vector<int> outGoingArcs;
	int size;
public:
	static const int SOURCE_NODE_ID;
	static const int SINK_NODE_ID;
	MPCNode(int id){
		this->id=id;
		this->lemonId=0;
		this->coverage=0;
		this->size=0;
	}

	double getCoverage() const {
		return coverage;
	}

	void setCoverage(double coverage) {
		this->coverage = coverage;
	}

	int getId() const {
		return id;
	}

	int getLemonId() const {
		return lemonId;
	}

	void setLemonId(int lemonId) {
		this->lemonId = lemonId;
	}

	void addIncomingArc(int arcId) {
		this->incomingArcs.push_back(arcId);
	}
	void addOutgoingArc(int arcId){
		this->outGoingArcs.push_back(arcId);
	}

	const vector<int>& getIncomingArcs() const {
		return incomingArcs;
	}

	const vector<int>& getOutGoingArcs() const {
		return outGoingArcs;
	}

	int getSize() const {
		return size;
	}

	void setSize(int size) {
		this->size = size;
	}
};
const int MPCNode::SOURCE_NODE_ID=-1;
const int MPCNode::SINK_NODE_ID=-2;

#endif /* MPCNODE_H_ */
