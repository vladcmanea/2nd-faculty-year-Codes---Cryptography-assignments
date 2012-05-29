#ifndef _timer_h
#define _timer_h

#include <windows.h>

class CPrecisionTimer
{
	LARGE_INTEGER _lFreq, _lStart;

public:
	CPrecisionTimer();
	void Start();
	double Stop();
};

#endif