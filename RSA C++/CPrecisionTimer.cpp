#include "CPrecisionTimer.h"

CPrecisionTimer::CPrecisionTimer()
{
	QueryPerformanceFrequency(&_lFreq);
}

void CPrecisionTimer::Start()
{
	QueryPerformanceCounter(&_lStart);
}

double CPrecisionTimer::Stop()
{
	// Return duration in seconds...
	LARGE_INTEGER lEnd;
	QueryPerformanceCounter(&lEnd);
	return (double(lEnd.QuadPart - _lStart.QuadPart) / _lFreq.QuadPart);
}
