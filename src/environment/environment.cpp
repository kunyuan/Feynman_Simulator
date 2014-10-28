//
//  environment.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/17/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "environment.h"
using namespace std;

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

Environment::Environment()
{
    Lat = nullptr;
}

Environment::~Environment()
{
    delete Lat;
}

bool Environment::BuildFromFile(string InputFile)
{
    ifstream ifs(InputFile, ios::in);
    ON_SCOPE_EXIT([&] {ifs.close(); });
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
    GET(_para, _L);

    Lat = new Lattice();

    GET(_para, Jcp);
    GET(_para, InitialBeta);
    GET(_para, DeltaBeta);
    GET(_para, Beta);
    GET(_para, Order);
    if (Order >= MAX_ORDER)
        ABORT("Order can not be bigger than " << MAX_ORDER);

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
    return true;
}

void Environment::SaveState()
{
}
