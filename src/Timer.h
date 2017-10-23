#ifndef TIMER_H
#define TIMER_H
/**
@file Timer.h
*/

#include <sys/time.h>
#include <iostream>

/**
 * @return returns current clock time in seconds 
 */
double tv_sec(struct timeval* tp){
	return ((tp->tv_sec + tp->tv_usec*1.0e-6));
}

/**
 * @class Timer
 * @brief times time in seconds
 */
class Timer
{
private:
	struct timeval start_tv;
	
public:
	/**
	 * @brief Starts timer
	 */
	void start(){
		gettimeofday(&start_tv, NULL);
	}
	/**
	 * @return time elapsed since last start
	 */
	double elapsed(){
		struct timeval elapsed_tv;
		gettimeofday(&elapsed_tv, NULL);
		elapsed_tv.tv_sec -= start_tv.tv_sec;
		elapsed_tv.tv_usec -= start_tv.tv_usec;
		return tv_sec(&elapsed_tv);
	}
};

#endif