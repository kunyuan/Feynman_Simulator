//
//  job.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/1/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "job.h"
#include <iostream>

using namespace std;

para::Job::Job(string inputfile)
{
    _Para.ParseFile(inputfile);
    GetPara(_Para, Type);
    if (TypeName.find(Type) == TypeName.end())
        ABORT("I don't know what is Job Type " << Type << "?");

    GetPara(_Para, DoesLoad);
    GetPara(_Para, PID);
    GetPara(_Para, WeightFile);
    GetPara(_Para, MessageFile);
    ParaFile = ToString(PID) + "_para.txt";
    StatisticsFile = ToString(PID) + "_statistics.npz";
    ConfigFile = ToString(PID) + "_config.txt";
    LogFile = ToString(PID) + ".log";
    InputFile = inputfile;
}