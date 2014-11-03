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
        Weight.BuildNew(weight::GW | weight::SigmaPolar, Para.Lat, Para.Beta, Para.Order);
    else {
        Weight.Load(_WeightFile(), weight::GW, Para.Lat, Para.Beta, Para.Order);
        Weight.BuildNew(weight::SigmaPolar, Para.Lat, Para.Beta, Para.Order);
    }
    Diag.BuildNew(Para.Lat, Weight.G, Weight.W);
    Grasshopper.BuildNew(Para, Diag, Weight);
    Scarecrow.BuildNew(Para, Diag, Weight);
    return true;
}

bool EnvMonteCarlo::Load()
{
    Para.Load(_ParaFile());
    Weight.Load(_WeightFile(), weight::GW | weight::SigmaPolar, Para.Lat, Para.Beta, Para.Order);
    Diag.Load(_ConfigFile(), Para.Lat, Weight.G, Weight.W);
    Grasshopper.BuildNew(Para, Diag, Weight);
    Scarecrow.Load(_StatisFile(), Para, Diag, Weight);
    return true;
}

void EnvMonteCarlo::Save()
{
    Para.Save(_ParaFile());
    Weight.Save(_WeightFile(), weight::GW | weight::SigmaPolar);
    Diag.Save(_ConfigFile());
    Scarecrow.Save(_StatisFile());
}

/**
*  Adjust everything according to new parameters, like new Beta, Jcp
*/
void EnvMonteCarlo::Reset(ParameterMap map)
{
    //    T = 1.0 / Beta;
    //    Sigma->Reset(Beta);
    //    Polar->Reset(Beta);
    //    G->Reset(Beta);
    //    W->Reset(Beta);
}
string EnvMonteCarlo::_ConfigFile()
{
    return ToString(PID) + "_config_env.txt";
}
