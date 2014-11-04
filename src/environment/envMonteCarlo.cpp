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
    _ParameterFile = ToString(PID) + "_para.txt";
    _GWweightFile = "GWweight";
    _WeightFile = ToString(PID) + "_statistics";
    _StatisticsFile = _WeightFile;
    _DiagramFile = ToString(PID) + "_diagram.txt";
}

bool EnvMonteCarlo::BuildNew(const std::string &InputFile, bool StartFromBare)
{
    //Read more stuff for the state of MC only
    Para.BuildNew(InputFile);
    if (StartFromBare)
        Weight.BuildNew(weight::GW | weight::SigmaPolar, Para);
    else {
        //Load GW weight from a global file shared by other MC processes
        Weight.Load(_GWweightFile, weight::GW, Para);
        Weight.BuildNew(weight::SigmaPolar, Para);
    }
    Diag.BuildNew(Para.Lat, Para.RNG, Weight.G, Weight.W);
    Grasshopper.BuildNew(Para, Diag, Weight);
    Scarecrow.BuildNew(Para, Diag, Weight);
    return true;
}
/**
*  Load() is used to continue an abrupt job, it loads GW weight and SigmaPolar weight from the same weight file
*
*  @return true if load succesfully
*/
bool EnvMonteCarlo::Load()
{
    Para.Load(_ParameterFile);
    Weight.Load(_StatisticsFile, weight::GW | weight::SigmaPolar, Para);
    Scarecrow.Load(_StatisticsFile, Para, Diag, Weight);
    Diag.Load(_DiagramFile, Para.Lat, Para.RNG, Weight.G, Weight.W);
    Grasshopper.BuildNew(Para, Diag, Weight);
    return true;
}

void EnvMonteCarlo::Save()
{
    Para.Save(_ParameterFile, "w");
    Weight.Save(_StatisticsFile, weight::GW | weight::SigmaPolar, "w");
    Scarecrow.Save(_StatisticsFile, "a"); // Save to the same file now
    Diag.Save(_DiagramFile, "w");
}

/**
*  Adjust everything according to new parameters, like new Beta, Jcp
*/
void EnvMonteCarlo::ReWeight(const Status &status)
{
    if (Para.Version >= status.Version)
        return;
    Para.Jcp = status.Jcp;
    Para.Beta = status.Beta;
    Para.T = 1.0 / Para.Beta;
    Weight.ReWeight(weight::GW | weight::SigmaPolar, Para);
    Grasshopper.ReWeight(Para);
    Scarecrow.ReWeight();
}
