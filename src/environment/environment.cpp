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
    Sigma = nullptr;
    Polar = nullptr;
    G = nullptr;
    W = nullptr;
}

Environment::~Environment()
{
    delete Sigma;
    delete Polar;
    delete G;
    delete W;
    delete Lat;
}

bool Environment::_ReadToPara(string file)
{
    ifstream ifs(file, ios::in);
    ON_SCOPE_EXIT([&] {ifs.close(); });
    if (!ifs.is_open())
        ABORT("Fail to open input file " << file);
    string temp;
    _para.clear();
    while (getline(ifs, temp))
        _para.push_back(temp);
    return true;
}

bool GetBool(const string source)
{
    if (source == "true" || source == "True")
        return true;
    else if (source == "false" || source == "False")
        return false;
    else
        ABORT("Fail to read DoesLoad");
}

bool Environment::BuildFromFile(string InputFile)
{
    if (!_ReadToPara(InputFile))
        return false;
    int type;
    GET(_para, type);
    Type = (JobType)type;

    string doesload;
    GET(_para, doesload);
    DoesLoad = GetBool(doesload);

    string startfrombare;
    GET(_para, startfrombare);
    StartFromBare = GetBool(startfrombare);

    GET(_para, PID);
    GET(_para, _L);
    GET(_para, Jcp);
    GET(_para, InitialBeta);
    GET(_para, DeltaBeta);
    GET(_para, FinalBeta);

    Beta = InitialBeta;
    T = 1.0 / Beta;

    GET(_para, Order);
    if (Order >= MAX_ORDER)
        ABORT("Order can not be bigger than " << MAX_ORDER);

    //make sure old Lat/Sigma/Polar/G/W are released before assign new memory
    delete Lat;
    Lat = new Lattice(_L);
    delete Sigma;
    Sigma = new Weight::Sigma(*Lat, Beta, Order);
    delete Polar;
    Polar = new Weight::Polar(*Lat, Beta, Order);
    delete G;
    G = new Weight::G(*Lat, Beta, Order);
    delete W;
    W = new Weight::W(*Lat, Beta, Order);
    return true;
}

string Environment::_StateFile()
{
    return ToString(PID) + "_state_env.txt";
}

string Environment::_GWweightFile()
{
    return "global_weight_env.npz";
}

string Environment::_LogFile()
{
    return ToString(PID) + ".log";
}

string Environment::_ControlFile()
{
    return "global_control.txt";
}

void Environment::SaveState(string Mode)
{
    auto mode = ios::out;
    if ((Mode) == "a")
        mode = ios::app;
    else if ((Mode) != "w")
        ABORT("I don't know what is the mode " << Mode << "?");

    ofstream ofs(_StateFile(), mode);
    if (!ofs.is_open())
        ABORT("Fail to open file " << _StateFile());
    ON_SCOPE_EXIT([&] {ofs.close(); });

    PUT(ofs, _L);
    PUT(ofs, Jcp);
    PUT(ofs, InitialBeta);
    PUT(ofs, DeltaBeta);
    PUT(ofs, FinalBeta);
    PUT(ofs, Beta);
    PUT(ofs, Order);
}

bool Environment::LoadState()
{
    if (!_ReadToPara(_StateFile()))
        return false;
    GET(_para, _L);
    Lat = new Lattice(_L);

    GET(_para, Jcp);
    GET(_para, InitialBeta);
    GET(_para, DeltaBeta);
    GET(_para, FinalBeta);
    GET(_para, Beta);
    T = 1.0 / Beta;
    GET(_para, Order);
    return true;
}

bool Environment::LoadGWweight()
{
    G->Load(_GWweightFile());
    W->Load(_GWweightFile());
    return true;
}

void Environment::SaveGWweight(string Mode)
{
    G->Save(_GWweightFile(), Mode);
    W->Save(_GWweightFile(), "a");
}
