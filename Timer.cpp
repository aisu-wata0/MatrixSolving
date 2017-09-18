#include <iostream>
#include <time.h>

namespace std{

double timeNow(void){
	struct timeval time;
	gettimeofday(&time, NULL);
	return ((double)(time.tv_sec*1000.0 + time.tv_usec/1000.0));
}

class Timer
{
private:
	double start_timer;
	
public:
	void start(){
		start_timer = timeNow();
	}
	
	double elapsed(){
		return total = timeNow() - start;
	}
};

} // namespace std