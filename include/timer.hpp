#ifndef timer_h
#define timer_h

class Timer {
public:
  void start();
  void stop();
  void get_elapsed_time();
  void get_elapsed_time(const char*);
private:
  double tStart;
  double tStop;
};

#endif /* timer_h */
