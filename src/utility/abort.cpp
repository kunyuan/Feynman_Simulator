//
//  abort.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 2/9/15.
//  Copyright (c) 2015 Kun Chen. All rights reserved.
//

#include "abort.h"
#include <csignal>

int InterruptHandler::__Signal = -1;

void InterruptHandler::__SignalHandler(int signum)
{
    LOG_INFO("Signal " << signum << " received. Interrupting program...");
    __Signal = signum;
    exit(signum);
}
void InterruptHandler::__DelayedSignalHandler(int signum)
{
    LOG_INFO("Delaying interrupt Signal " << signum);
    __Signal = signum;
}

InterruptHandler::InterruptHandler()
{
    __IsDelaying = false;
    __Signal = -1;
    signal(SIGINT, __SignalHandler);
    signal(SIGTERM, __SignalHandler);
}

InterruptHandler::~InterruptHandler()
{
    signal(SIGINT, SIG_DFL);
    signal(SIGTERM, SIG_DFL);
}
void InterruptHandler::Delay()
{
    signal(SIGINT, __DelayedSignalHandler);
    signal(SIGTERM, __DelayedSignalHandler);
    __IsDelaying = true;
}

void InterruptHandler::Resume()
{
    signal(SIGINT, __SignalHandler);
    signal(SIGTERM, __SignalHandler);
    if (__IsDelaying && (__Signal == SIGINT || __Signal == SIGTERM)) {
        __SignalHandler(__Signal);
        __Signal = -1;
        __IsDelaying = false;
    }
}
