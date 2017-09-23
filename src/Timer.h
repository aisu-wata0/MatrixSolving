#ifndef TIMER_H
#define TIMER_H

#include <sys/time.h>

/**
 * @return returns current clock time in seconds 
 */
double timestamp(void){
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return ((double)(tp.tv_sec + tp.tv_usec/1000000.0));
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
		start_timer = timestamp();
	}
	/**
	 * @return time elapsed since last start
	 */
	double elapsed(){
		return timestamp() - start_timer;
	}
};

#endif // TIMER_H