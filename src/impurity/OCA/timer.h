// @Copyright 2007 Kristjan Haule
// 
#ifndef trtimer
#define trtimer
#include <cstdio>
#include <cstdlib>
#include <sys/types.h>
#include <sys/times.h>
#include <unistd.h> 
#include <ctime>

class Timer{
public:
  void start() {times(&_start);}
  void stop() {times(&_end);}
  char* ToSeconds();
  char* ToTicks();
  double toSeconds();
  Timer() : CLK_TCK_(sysconf(_SC_CLK_TCK)){};
private:
  tms _start, _end;
  char bufferT[100], bufferS[100];
  long CLK_TCK_;
};

char* Timer::ToSeconds()
{
  double utime = double(_end.tms_utime - _start.tms_utime) / CLK_TCK_;
  double stime = double(_end.tms_stime - _start.tms_stime) / CLK_TCK_;
  sprintf(bufferS, "%.2fu %.2fs", utime, stime);
  return bufferS;
}

double Timer::toSeconds()
{
  double utime = double(_end.tms_utime - _start.tms_utime) / CLK_TCK_;
  double stime = double(_end.tms_stime - _start.tms_stime) / CLK_TCK_;
  return utime+stime;
}

char* Timer::ToTicks()
{
  clock_t utime = _end.tms_utime - _start.tms_utime;
  clock_t stime = _end.tms_stime - _start.tms_stime;
  sprintf(bufferT, "%du %ds", static_cast<int>(utime), static_cast<int>(stime));
  return bufferT;
}
#endif
