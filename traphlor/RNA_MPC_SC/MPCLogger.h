/*
 * MPCLogger.h
 *
 *  Created on: Feb 11, 2014
 *      Author: ahmedsobih
 */

#ifndef MPCLOGGER_H_
#define MPCLOGGER_H_

#include "MPCHeaders.h"
class MPCLogger {
private:
	static const bool debug = false;
	ofstream logFile;
	static string getTime(){
		struct tm *aTime ;
		time_t theTime = time(NULL);
		aTime = localtime(&theTime);
		char timeBuffer[80];
		strftime(timeBuffer,80,"%G%m%d%H%M%S",aTime);
		return timeBuffer;
	}

public:
	static void log(string msg){
		if(debug)
			cout<< msg<<endl;
	}
	MPCLogger();
	virtual ~MPCLogger();
};

#endif /* MPCLOGGER_H_ */
