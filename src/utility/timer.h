//
// Created by Kun Chen on 10/5/14.
// Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __timer_H_
#define __timer_H_

#include <ctime>
#include <iosfwd>
#include <iomanip>

class timer {
    friend std::ostream& operator<<(std::ostream& os, timer& t);

private:
    bool running;
    clock_t start_clock;
    time_t start_time;
    double acc_time;

    double elapsed_time();

public:
    // 'running' is initially false.  A timer needs to be explicitly started
    // using 'start' or 'restart'
    timer()
        : running(false)
        , start_clock(0)
        , start_time(0)
        , acc_time(0)
    {
    }

    void start(const char* msg = 0);
    void restart(const char* msg = 0);
    void stop(const char* msg = 0);
    void check(const char* msg = 0);
    /**
    *  check whether Interval seconds has passed since the timer is started. Once it is true, the timer will be automatically restarted
    *
    *  @param Interval time interval in Seconds
    *
    */
    bool check(time_t Interval);

}; // class timer

int TestTimer();

#endif //__timer_H_
