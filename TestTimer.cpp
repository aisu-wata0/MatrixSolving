#include <iostream>
#include <boost/date_time/posix_time/ptime.hpp>
#include <boost/date_time/microsec_time_clock.hpp>

using namespace std;
using namespace boost;

class TestTimer
{
private:
	string name;
	posix_time::ptime start;
	
public:
	TestTimer(const std::string & name) : name(name),
		start(date_time::microsec_clock<posix_time::ptime>::local_time())
	{
		cout <<"\n>>>== "<< name <<" started"<< endl;
	}

	~TestTimer()
	{

		posix_time::ptime now(date_time::microsec_clock<posix_time::ptime>::local_time());
		posix_time::time_duration d = now - start;
		
		cout <<"\n<<<== "<< name <<" completed in " << d.total_milliseconds() / 1000.0 <<
			" seconds" << endl;
	}
};