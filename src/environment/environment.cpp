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
}

string Environment::_ParaFile()
{
    return ToString(PID) + "_para_env.txt";
}

string Environment::_WeightFile()
{
    return ToString(PID) + "_weight_env.npz";
}

string Environment::_ControlFile()
{
    return "global_control.txt";
}

string Environment::_StatisFile()
{
    return ToString(PID) + "_statis_env.npz";
}
