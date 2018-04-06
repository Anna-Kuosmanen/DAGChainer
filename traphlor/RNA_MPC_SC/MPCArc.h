/*
 * MPCArc.h
 *
 *  Created on: Feb 11, 2014
 *      Author: ahmedsobih
 */

#ifndef MPCARC_H_
#define MPCARC_H_

#define int64_t_MAX std::numeric_limits<int64_t>::max()

class MPCArc {
private:
	int id;
	int source;
	int target;
	int lemonId;
	double cost;
	double coverage;
	int lowerBound;
	long upperBound;
	int subpathId;
public:
	MPCArc(int id, int source, int destination){
		this->id=id;
		this->source=source;
		this->target=destination;
		this->lemonId=0;
		this->cost=0;
		this->coverage=0;
		this->lowerBound=0;
		this->upperBound=int64_t_MAX;
		subpathId=-1;
	}

	int getTarget() const {
		return target;
	}

	void setTarget(int target) {
		this->target = target;
	}

	int getId() const {
		return id;
	}

	void setId(int id) {
		this->id = id;
	}

	int getLemonId() const {
		return lemonId;
	}

	void setLemonId(int lemonId) {
		this->lemonId = lemonId;
	}

	int getSource() const {
		return source;
	}

	void setSource(int source) {
		this->source = source;
	}

	double getCost() const {
		return cost;
	}

	void setCost(double cost) {
		this->cost = cost;
	}

	double getCoverage() const {
		return coverage;
	}

	void setCoverage(double coverage) {
		this->coverage = coverage;
	}

	int getLowerBound() const {
		return lowerBound;
	}

	void increaseLowerBound() {
		this->lowerBound +=1;
	}

	long getUpperBound() const {
		return upperBound;
	}

	void setUpperBound(long value) {
		this->upperBound = value;
	}

	int getSubpathId() const {
		return subpathId;
	}

	void setSubpathId(int subpathId) {
		this->subpathId = subpathId;
	}
};

#endif /* MPCARC_H_ */
