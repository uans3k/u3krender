#include "stdafx.h"
#include "Timer.h"


Timer::Timer()
{
	fps = 30;
	frameCount = 1000 / fps;
}

Timer::Timer(int fps)
{
	this->fps = fps;
	frameCount = 1000 / fps;
}


Timer::~Timer()
{
}

int Timer::start()
{
	startCount = GetTickCount();
	return startCount;
}

int Timer::wait()
{	
	while ((GetTickCount() - startCount) < frameCount);
	return GetTickCount();
}


