//
// Created by Kun Chen on 10/5/14.
// Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "timer.h"
#include "logger.h"

//===========================================================================
// Return the total time that the timer has been in the "running"
// state since it was first "started" or last "restarted".  For
// "short" time periods (less than an hour), the actual cpu time
// used is reported instead of the elapsed time.

double timer::elapsed_time()
{
    time_t acc_sec = time(0) - start_time;
    if (acc_sec < 3600)
        return (clock() - start_clock) / (1.0 * CLOCKS_PER_SEC);
    else
        return (1.0 * acc_sec);

} // timer::elapsed_time

//===========================================================================
// Start a timer.  If it is already running, let it continue running.
// Print an optional message.

void timer::start(const char* msg)
{
    // Print an optional message, something like "Starting timer t";
    if (msg)
        LOG_INFO(msg << std::endl);

    // Return immediately if the timer is already running
    if (running)
        return;

    // Set timer status to running and set the start time
    running = true;
    start_clock = clock();
    start_time = time(0);

} // timer::start

//===========================================================================
// Turn the timer off and start it again from 0.  Print an optional message.

void timer::restart(const char* msg)
{
    // Print an optional message, something like "Restarting timer t";
    if (msg)
        LOG_INFO(msg << std::endl);

    // Set timer status to running, reset accumulated time, and set start time
    running = true;
    acc_time = 0;
    start_clock = clock();
    start_time = time(0);

} // timer::restart

//===========================================================================
// Stop the timer and print an optional message.

void timer::stop(const char* msg)
{
    // Print an optional message, something like "Stopping timer t";
    if (msg)
        LOG_INFO(msg << std::endl);

    // Compute accumulated running time and set timer status to not running
    if (running)
        acc_time += elapsed_time();
    running = false;

} // timer::stop

//===========================================================================
// Print out an optional message followed by the current timer timing.

void timer::check(const char* msg)
{
    // Print an optional message, something like "Checking timer t";
    if (msg)
        LOG_INFO(msg << " : ");

    LOG_INFO("Elapsed time [" << std::setiosflags(std::ios::fixed)
                              << std::setprecision(2)
                              << acc_time + (running ? elapsed_time() : 0) << "] seconds");

} // timer::check

bool timer::check(time_t Interval)
{
    time_t acc_sec = time(0) - start_time;
    if (acc_sec < Interval)
        return false;
    else {
        restart();
        return true;
    }
}

//===========================================================================
// Allow timers to be printed to ostreams using the syntax 'os << t'
// for an ostream 'os' and a timer 't'.  For example, "cout << t" will
// print out the total amount of time 't' has been "running".

std::ostream& operator<<(std::ostream& os, timer& t)
{
    os << std::setprecision(2) << std::setiosflags(std::ios::fixed)
       << t.acc_time + (t.running ? t.elapsed_time() : 0);
    return os;
}

//===========================================================================
