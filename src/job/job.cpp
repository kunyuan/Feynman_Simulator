//
//  parameter.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/3/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "job.h"
#include "../utility/abort.h"
#define GET(is, thing)                               \
    {                                                \
        string temp;                                 \
        (*is) >> thing;                              \
        getline((*is), temp);                        \
        if (!(*is).good())                           \
            ABORT("Fail to read " << #thing << "!"); \
    }

using namespace std;
using namespace Jobs;

void OpenFile(string InputFile, jobstream ifs)
{
    (*ifs).open(InputFile, ios::in);
    if (!(*ifs).is_open())
        ABORT("Fail to open input file " << InputFile);
}

JobType GetJobsType(std::string InputFile)
{
    jobstream ifs;
    OpenFile(InputFile, ifs);
    int type;
    GET(ifs, type);
    return (JobType)type;
}

JobsBase::JobsBase(string InputFile)
{
    jobstream ifs;
    OpenFile(InputFile, ifs);
    int type;
    GET(ifs, type);
    Type = (JobType)type;
    GET(ifs, PID);
    GET(ifs, L);
    GET(ifs, Jcp);
    GET(ifs, InitialBeta);
    GET(ifs, DeltaBeta);
    GET(ifs, Beta);
    GET(ifs, Order);
    string doesload;
    GET(ifs, doesload);
    if (doesload == "true") {
        GET(ifs, StateFile);
        DoesLoad = true;
    }
    else if (doesload == "false") {
        DoesLoad = false;
    }
    else
        ABORT("Fail to read DoesLoad");
}

void JobsBase::TestJob()
{
    Beta = 1.0;
    Order = 1;
    Type = MC;
}

JobsMC::JobsMC(string InputFile)
    : JobsBase(InputFile)
{
}
