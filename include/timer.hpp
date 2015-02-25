#ifndef timer_h
#define timer_h

#include <unistd.h>
#include <sys/time.h>
#include <iostream>

// Return time in seconds since the Unix epoch
double timer(void);

class Timer {
public:
  void start();
  void stop();
  void get_elapsed_time();
private:
  double tStart;
  double tStop;
};

inline void Timer::start() {
  tStart = timer();
}

inline void Timer::stop() {
  tStop = timer();
}

inline void Timer::get_elapsed_time() {
  std::cout << "Elapsed time : " << tStop-tStart << std::endl;
}

inline double timer(void)
{
  double time;
  struct timeval tv;
  gettimeofday(&tv, NULL);
  time = (double)tv.tv_sec + (double)tv.tv_usec/1.e6;
  return time;
}


#endif /* timer_h */
