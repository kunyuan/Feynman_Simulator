//
//  job.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/1/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "job.h"
#include "utility/dictionary.h"
#include <iostream>

using namespace std;

para::Job::Job(string inputfile)
{

    Dictionary _Para;
    _Para.Load(inputfile);
    _Para = _Para.Get<Dictionary>("Job");
    GET(_Para, Type);
    if (TypeName.find(Type) == TypeName.end())
        ABORT("I don't know what is Job Type " << Type << "?");

    GET(_Para, DoesLoad);
    GET(_Para, Sample);
    GET(_Para, PID);
    GET(_Para, Sample);
    GET(_Para, WeightFile);
    GET(_Para, MessageFile);
    string Prefix = ToString(PID) + "_" + string(Type);
    ParaFile = Prefix + "_para";
    StatisticsFile = Prefix + "_statis";
    LogFile = Prefix + ".log";
    InputFile = inputfile;
}