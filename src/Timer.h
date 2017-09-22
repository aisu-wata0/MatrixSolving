#ifndef TIMER_H
#define TIMER_H

#include <sys/time.h>

double gettime(void){
	timeval time;
	gettimeofday(&time, NULL);
	return ((double)(time.tv_sec + time.tv_usec/1000000.0));
}

/**
 * @class Timer
 * @brief times time in seconds
 */
class Timer
{
private:
	double start_timer;
	
public:
	/**
	 * @brief Starts timer
	 */
	void start(){
		start_timer = gettime();
	}
	/**
	 * @return time elapsed since last start
	 */
	double elapsed(){
		return gettime() - start_timer;
	}
};

#endif // TIMER_H