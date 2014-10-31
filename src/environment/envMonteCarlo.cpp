//
//  envMonteCarlo.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/18/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "environment.h"

using namespace std;
EnvMonteCarlo::EnvMonteCarlo()
{
    Version = 0;
}

EnvMonteCarlo::~EnvMonteCarlo()
{
}

bool EnvMonteCarlo::BuildFromFile(string InputFile)
{
    Environment::BuildFromFile(InputFile);
    LOGGER_CONF(_LogFile(), "MC", Logger::file_on | Logger::screen_on, INFO, INFO);

    if (Type != MC)
        ABORT("The job input file is for Monte Carlo!");

    //Read more parameters
    GET(_para, Toss);
    GET(_para, Sample);
    GET(_para, Sweep);
    GET(_para, Seed);
    GET(_para, WormSpaceReweight);
    _ReadOrderWeight();

    //Initialize utilies for MC simulations
    Counter = 0;
    WormWeight = Weight::Worm();

    Diag.Build(Lat, G, W);

    //Initialize random number generator
    this->RNG.Reset(Seed);

    //Add estimators
    cEstimator.AddEstimator("1");
    cEstimator.AddEstimator("2");
    rEstimator.AddEstimator("2");
    rEstimator.AddEstimator("3");

    Load(); //call load to initialize staff in EnvMonteCarlo
    return true;
}

/**
*  Adjust everything according to new parameters, like new Beta, Jcp
*/
void EnvMonteCarlo::Reset()
{
    T = 1.0 / Beta;
    Sigma->Reset(Beta);
    Polar->Reset(Beta);
    G->Reset(Beta);
    W->Reset(Beta);
}

bool EnvMonteCarlo::Load()
{
    if (StartFromBare) {
        //TODO: initialize GW with bare weight
    }
    else
        LoadGWweight();

    if (DoesLoad) {
        LoadState();
        LoadConfig();
        LoadStatis();
    }
    return true;
}

void EnvMonteCarlo::Save()
{
    SaveState("w");
    SaveConfig("w");
    SaveStatis("w");
}

void EnvMonteCarlo::SaveState(string Mode)
{
    //Save parameters to .state.env file
    Environment::SaveState(Mode);
    ofstream ofs(_StateFile(), ios::app);
    //append to the state file which has already been created by Environment::SaveState()
    ON_SCOPE_EXIT([&] {ofs.close(); });
    if (!ofs.is_open())
        ABORT("Fail to open file " << _StateFile());
    PUT(ofs, Counter);
    PUT(ofs, Toss);
    PUT(ofs, Sample);
    PUT(ofs, Sweep);
    PUT(ofs, WormSpaceReweight);
    _WriteOrderWeight(ofs);
    PUT(ofs, RNG);
}

bool EnvMonteCarlo::LoadState()
{
    Environment::LoadState();
    GET(_para, Counter);
    GET(_para, Toss);
    GET(_para, Sample);
    GET(_para, Sweep);
    GET(_para, WormSpaceReweight);
    _ReadOrderWeight();
    GET(_para, this->RNG);
    return true;
}

void EnvMonteCarlo::SaveStatis(string Mode)
{
    //Save statistics to .statis.env
    cEstimator.SaveStatistics(_StatisFile(), Mode);
    rEstimator.SaveStatistics(_StatisFile(), "a");

    Sigma->Save(_StatisFile(), "a");
    Polar->Save(_StatisFile(), "a");
}

bool EnvMonteCarlo::LoadStatis()
{
    cEstimator.LoadStatistics(_StatisFile());
    rEstimator.LoadStatistics(_StatisFile());

    Sigma->Load(_StatisFile());
    Polar->Load(_StatisFile());
    return true;
}

void EnvMonteCarlo::SaveConfig(string Mode)
{
    //Save diagram to .config.env file
    Diag.SaveConfig(_ConfigFile(), Mode);
}

bool EnvMonteCarlo::LoadConfig()
{
    //Save diagram to .config.env file
    Diag.LoadConfig(_ConfigFile());
    return true;
}

string EnvMonteCarlo::_ConfigFile()
{
    return ToString(PID) + "_config_env.txt";
}

string EnvMonteCarlo::_StatisFile()
{
    return ToString(PID) + "_statis_env.npz";
}

void EnvMonteCarlo::_WriteOrderWeight(ostream &os, char sep)
{
    os << OrderWeight[0];
    for (int i = 1; i < Order; i++) {
        os << sep << OrderWeight[i];
    }
    os << "   #OrderWeight" << endl;
    if (os.fail())
        ABORT("Read OrderWeight fails!");
}

void EnvMonteCarlo::_ReadOrderWeight(char sep)
{
    stringstream ss(_para.front());
    char sepchar;
    ss >> OrderWeight[0];
    for (int i = 1; i < Order; i++) {
        ss >> sepchar;
        if (sepchar != sep)
            ABORT("Sep char " << sepchar << " is not expected as the separator. I will expect " << sep);
        ss >> OrderWeight[i];
    }
    if (ss.fail())
        ABORT("Read OrderWeight fails!");
    _para.pop_front();
}
