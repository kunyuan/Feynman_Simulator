//
//  enDyson.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/18/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "environment.h"

EnvDyson::EnvDyson(int pid)
    : Environment(pid)
{
    _ParameterFile = ToString(PID) + "_para.txt";
    _GWweightFile = "GWweight.npz";
    _WeightFile = ToString(PID) + "_statistics.npz";
    _StatisticsFile = _WeightFile;
}

bool EnvDyson::BuildNew(const string &InputFile, bool StartFromBare)
{
    Para.BuildNew(InputFile);
    if (StartFromBare) {
        LOG_INFO("Build Dyson parameters from " << _ParameterFile);
        Weight.BuildNew(weight::GW | weight::SigmaPolar, Para);
    }
    else {
        //Load GW weight from a global file shared by other MC processes
        Weight.Load(_GWweightFile, weight::GW, Para);
        Weight.Load(_WeightFile, weight::SigmaPolar, Para);
    }
    return true;
}

bool EnvDyson::CanBeLoad()
{
    return DoesFileExist(_ParameterFile);
}

bool EnvDyson::Load()
{
    LOG_INFO("Loading Dyson environment from " << _ParameterFile);
    Para.Load(_ParameterFile);
    Weight.Load(_GWweightFile, weight::GW, Para);
    Weight.Load(_WeightFile, weight::SigmaPolar, Para);
    return true;
}

void EnvDyson::Save()
{
    LOG_INFO("Saving Dyson environment to " << _ParameterFile);
    Para.Save(_ParameterFile, "w");
    Weight.Save(_GWweightFile, weight::GW, "w");
    Weight.Save(_WeightFile, weight::SigmaPolar, "w");
    Para.Version++;
    LOG_INFO("Dyson Version becomes " << Para.Version);
    Para.GetStatus().Save();
}
