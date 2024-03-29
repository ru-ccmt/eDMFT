// @Copyright 2007 Kristjan Haule
// 
#ifndef __TIMER__
#define __TIMER__
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <sys/types.h>
#include <sys/times.h>
#include <unistd.h>

/*
class Timer {
  timespec starttime;
  double elapsed_time;		// elapsed time in seconds

public:
  Timer() : elapsed_time(0) {}

  void reset() { elapsed_time = 0; }
  void start() { clock_gettime(CLOCK_MONOTONIC, &starttime); }
  void stop() {
    timespec endtime;
    clock_gettime(CLOCK_MONOTONIC, &endtime);
    elapsed_time += (double)(endtime.tv_sec - starttime.tv_sec)
      + ((double)(endtime.tv_nsec - starttime.tv_nsec))/(1.0e9);
  }
  double elapsed() const { return elapsed_time; }
};
*/

#ifdef _TIME
class Timer{
  double elapsed_time;          // elapsed time in seconds
public:
  Timer() : CLK_TCK_(sysconf(_SC_CLK_TCK))
  { elapsed_time=0; }
  void reset() {elapsed_time=0;}
  void start() {times(&_start);}
  void stop() {
    times(&_end);
    double utime = double(_end.tms_utime - _start.tms_utime) / CLK_TCK_;
    double stime = double(_end.tms_stime - _start.tms_stime) / CLK_TCK_;
    elapsed_time +=  utime+stime;
  }
  double elapsed() const { return elapsed_time; }
private:
  tms _start, _end;
  long CLK_TCK_;
};
#else
class Timer{
public:
  Timer() {};
  void reset() {};
  void start() {};
  void stop() {};
  double elapsed() {return 0.0;}
};
#endif

#endif /* __TIMER__ */
