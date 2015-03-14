#include "timer.hpp"

#include <unistd.h>
#include <sys/time.h>
#include <iostream>

// Return time in seconds since the Unix epoch
static double timer() {
double time;
  struct timeval tv;
  gettimeofday(&tv, NULL);
  time = (double)tv.tv_sec + (double)tv.tv_usec/1.e6;
  return time;
}

void Timer::start() {
  tStart = timer();
}

void Timer::stop() {
  tStop = timer();
}

void Timer::get_elapsed_time() {
std::cout << "Elapsed time : " << tStop-tStart << std::endl;
}

void Timer::get_elapsed_time(const char* msg) {
std::cout << msg << " : "
<< tStop-tStart << " seconds" << std::endl;
}

