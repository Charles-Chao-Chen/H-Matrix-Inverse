#ifndef timer_h
#define timer_h

class Timer {
public:
  void start();
  void stop();
  void show_elapsed_time();
  void show_elapsed_time(const char*);
  double get_elapsed_time();
private:
  double tStart;
  double tStop;
};

#endif /* timer_h */
