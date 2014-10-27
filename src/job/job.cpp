//
//  parameter.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/3/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "job.h"
#include "../utility/abort.h"
using namespace std;

#define GET(para, thing)                             \
    {                                                \
        stringstream ss(para.front());               \
        (ss) >> thing;                               \
        if (!(ss).good())                            \
            ABORT("Fail to read " << #thing << "!"); \
        para.pop_front();                            \
    }

JobType GetJobsType(std::string InputFile)
{
    ifstream ifs(InputFile, ios::in);
    if (!ifs.is_open())
        ABORT("Fail to open input file " << InputFile);
    int type;
    ifs >> type;
    if (!ifs.good())
        ABORT("Fail to read JobType!");
    ifs.close();
    return (JobType)type;
}

JobsBase::JobsBase(string InputFile)
{
    ifstream ifs(InputFile, ios::in);
    if (!ifs.is_open())
        ABORT("Fail to open input file " << InputFile);
    string temp;
    _para.clear();
    while (getline(ifs, temp))
        _para.push_back(temp);

    int type;
    GET(_para, type);
    Type = (JobType)type;
    GET(_para, PID);
    GET(_para, L);
    GET(_para, Jcp);
    GET(_para, InitialBeta);
    GET(_para, DeltaBeta);
    GET(_para, Beta);
    GET(_para, Order);
    string doesload;
    GET(_para, doesload);
    if (doesload == "true") {
        GET(_para, StateFile);
        DoesLoad = true;
    }
    else if (doesload == "false") {
        StateFile = "";
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
    GET(_para, Toss);
    GET(_para, Sample);
    GET(_para, Sweep);
    GET(_para, Seed);
    GET(_para, WormSpaceReweight);
}

JobsMC::JobsMC()
{
    TestJob();
}

JobsDyson::JobsDyson(string InputFile)
    : JobsBase(InputFile)
{
}

JobsDyson::JobsDyson()
{
    TestJob();
}