//
//  envMonteCarlo.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/18/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "environment.h"

using namespace std;
EnvMonteCarlo::EnvMonteCarlo(int pid)
    : Environment(pid)
{
}

bool EnvMonteCarlo::BuildNew(const string &InputFile, bool StarFromBare)
{
    //Read more stuff for the state of MC only
    Para.BuildNew(InputFile);
    if (StarFromBare)
        Weight.BuildNew(weight::GW | weight::SigmaPolar, Para);
    else {
        Weight.Load(_WeightFile(), weight::GW, Para);
        Weight.BuildNew(weight::SigmaPolar, Para);
    }
    Diag.BuildNew(Para.Lat, Weight.G, Weight.W);
    Grasshopper.BuildNew(Para, Diag, Weight);
    Scarecrow.BuildNew(Para, Diag, Weight);
    return true;
}

bool EnvMonteCarlo::Load()
{
    Para.Load(_ParaFile());
    Weight.Load(_WeightFile(), weight::GW | weight::SigmaPolar, Para);
    Diag.Load(_ConfigFile(), Para.Lat, Weight.G, Weight.W);
    Grasshopper.BuildNew(Para, Diag, Weight);
    Scarecrow.Load(_StatisFile(), Para, Diag, Weight);
    return true;
}

void EnvMonteCarlo::Save()
{
    Para.Save(_ParaFile(), "w");
    Weight.Save(_WeightFile(), weight::GW | weight::SigmaPolar, "w");
    Diag.Save(_ConfigFile(), "w");
    Scarecrow.Save(_StatisFile(), "w");
}

/**
*  Adjust everything according to new parameters, like new Beta, Jcp
*/
void EnvMonteCarlo::ReWeight(const State &state)
{
    if (Para.Version >= state.Version)
        return;
    Para.Jcp = state.Jcp;
    Para.Beta = state.Beta;
    Para.T = 1.0 / Para.Beta;
    Weight.ReWeight(weight::GW | weight::SigmaPolar, Para);
    Grasshopper.ReWeight(Para);
    Scarecrow.ReWeight();
}

string EnvMonteCarlo::_ConfigFile()
{
    return ToString(PID) + "_config_env.txt";
}
