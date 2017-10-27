#ifndef TIMER_H
#define TIMER_H
/**
@file Timer.h
*/

#include <sys/time.h>
#include <iostream>
#include <cstddef>

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


template<size_t size>
class Chronometer
{
public:
	struct timeval clock[size];
	clock_t n_clock[size];
	size_t c;
	
	Chronometer(){
		init();
	}
	
	void init(){
		gettimeofday(&clock[0], NULL);
		n_clock[0] = std::clock();
		c = 0;
	}
	
	double tick(double* cl){
		struct timeval elapsed_tv;
		c = (c + 1) % (ptrdiff_t)size;
		gettimeofday(&clock[c], NULL);
		elapsed_tv.tv_sec = clock[c].tv_sec - clock[c-1].tv_sec;
		elapsed_tv.tv_usec = clock[c].tv_usec - clock[c-1].tv_usec;
		
		n_clock[c] = std::clock();
		*cl = (double)((n_clock[c] - n_clock[c-1]) / CLOCKS_PER_SEC);
		return tv_sec(&elapsed_tv);
	}
	
	void tick(){
		double dummy;
		tick(&dummy);
	}
};

#endif