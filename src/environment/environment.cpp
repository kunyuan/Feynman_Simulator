//
//  environment.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/17/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "environment.h"
using namespace std;

Environment::Environment(int pid)
    : PID(pid)
{
    _ParameterFile = ToString(PID) + "_para.txt";
    _GWweightFile = "GWweight.npz";
    _WeightFile = ToString(PID) + "_statistics.npz";
    _StatisticsFile = _WeightFile;
}