#pragma once
class Timer
{
public:
	Timer();
	Timer(int fps);
	~Timer();
	int start();
	int wait();
private:
	long startCount;
	long flowCount;
	long frameCount;
	int fps;
};

