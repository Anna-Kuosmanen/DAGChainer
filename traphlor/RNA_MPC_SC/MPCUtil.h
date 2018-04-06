/*
 * MPCUtil.h
 *
 *  Created on: Feb 11, 2014
 *      Author: ahmedsobih
 */

#include "MPCHeaders.h"

using namespace std;


#ifndef MPCUTIL_H_
#define MPCUTIL_H_
class MPCUtil {

public:
	static const string SPACE_SEPARATOR;
	static const string TAB_SEPARATOR;
	MPCUtil();
	virtual ~MPCUtil();

	static int getIntValue(string str){
		return atoi(str.c_str());
	}
	static double getDoubleValue(string str){
		return atof(str.c_str());
	}
	static string getStringValue(int num){
		return static_cast<ostringstream*>( &(ostringstream() << num) )->str();
	}
	static string getStringValue(long num){
		return static_cast<ostringstream*>( &(ostringstream() << num) )->str();
	}static string getStringValue(double num){
		return static_cast<ostringstream*>( &(ostringstream() << num) )->str();
	}
	static bool isIncluded(string a, string b){
		if(a=="" || b=="")
			return false;
		if(a.find(b)!=string::npos)
			return true;
		if(b.find(a)!=string::npos)
			return true;
		return false;
	}
	static double calcLogCost(double a, double b){
//		double sum=a+b;
//		return -log(1-abs((a-b)/sum));
		double cost = abs(a-b);
		if(cost == 0)
			cost = 1;
		return cost;
	}
	static string merge(string a, string b){
		string merged="";
		//cout<<"A is: "<<a<<"  B is: "<<b<<endl;
		vector<string> splittedB;
		vector<string> splittedA;
		split(splittedB, b, is_any_of(MPCUtil::SPACE_SEPARATOR));
		split(splittedA, a, is_any_of(MPCUtil::SPACE_SEPARATOR));
		if(splittedA[0]==splittedB[0] || splittedA[splittedA.size()-1]==splittedB[splittedB.size()-1])
			return merged;
		if(splittedA[splittedA.size()-1]==splittedB[0]){
			merged=a+b.substr(splittedB[0].size(),b.size()-1);
			//cout<<"Merged Path= "<<merged<<endl;
			return merged;
		}
		for(unsigned i=0;i<splittedA.size();i++){
			if(splittedA[i]==splittedB[0]){
				if(splittedB.size()<=splittedA.size()-i){
					return "";
				}else{
					int counterB=0;
					for(unsigned counterA=i;counterA<splittedA.size();counterA++,counterB++){
						if(splittedA[counterA]!=splittedB[counterB])
							return "";
					}
					merged=a+" "+b.substr(b.find(splittedB[counterB]), b.size()-1);
					//cout<<"Merged Path= "<<merged<<endl;
					return merged;
				}
			}
		}
		return merged;
	}
	static string mergePath(string a, string b){
		string merged=merge(a,b);
		if(merged=="")
			merge(b,a);
		return merged;
	}
};
const string MPCUtil::SPACE_SEPARATOR=" ";
const string MPCUtil::TAB_SEPARATOR="\t";
#endif /* MPCUTIL_H_ */
