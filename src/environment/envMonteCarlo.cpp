//
//  envMonteCarlo.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/18/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "environment.h"

using namespace std;
using namespace para;

EnvMonteCarlo::EnvMonteCarlo(const para::Job& job, bool IsAllTauSymmetric)
    : Job(job)
    , Weight(IsAllTauSymmetric)
{
}

bool EnvMonteCarlo::BuildNew()
{
    LOGGER_CONF(Job.LogFile, Job.Type, Logger::file_on | Logger::screen_on, INFO, INFO);
    //Read more stuff for the state of MC only
    Para.BuildNew(Job.InputFile);
    Para.Save(Job.InputFile, "w"); //save a copy of new para file
    //Load GW weight from a global file shared by other MC processes
    Weight.Load(Job.WeightFile, weight::GW, Para);
    Weight.BuildNew(weight::SigmaPolar, Para);
    Diag.BuildNew(Para.Lat, *Weight.G, *Weight.W);
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
    LOGGER_CONF(Job.LogFile, Job.Type, Logger::file_on | Logger::screen_on, INFO, INFO);
    Para.Load(Job.ParaFile);
    Weight.Load(Job.StatisticsFile, weight::GW | weight::SigmaPolar, Para);
    Diag.Load(Job.ConfigFile, Para.Lat, *Weight.G, *Weight.W);
    Grasshopper.BuildNew(Para, Diag, Weight);
    Scarecrow.Load(Job.StatisticsFile, Para, Diag, Weight);
    return true;
}

void EnvMonteCarlo::Save()
{
    Para.Save(Job.ParaFile, "w");
    Weight.Save(Job.StatisticsFile, weight::GW | weight::SigmaPolar, "w");
    Diag.Save(Job.ConfigFile, "w");
    Scarecrow.Save(Job.StatisticsFile, "a"); // Save to the same statis file as weight
}
void EnvMonteCarlo::DeleteSavedFiles()
{
    system(("rm " + Job.ParaFile).c_str());
    system(("rm " + Job.StatisticsFile).c_str());
    system(("rm " + Job.WeightFile).c_str());
    system(("rm " + Job.ConfigFile).c_str());
}

/**
*  Adjust everything according to new parameters, like new Beta, Jcp
*/
bool EnvMonteCarlo::ListenToMessage()
{
    LOG_INFO("Start reweighting...");
    Message Message_;
    if (!Message_.Load())
        return false;
    if (Para.Version >= Message_.Version) {
        LOG_INFO("Status has not been updated yet since the last reweighting!");
        return false;
    }
    Para.UpdateWithMessage(Message_);
    Weight.Load(Job.WeightFile, weight::GW, Para);
    Weight.ReWeight(weight::GW | weight::SigmaPolar, Para);
    Grasshopper.ReWeight(Para);
    Scarecrow.ReWeight();
    LOG_INFO("Reweighted to:\n" << Message_.PrettyString());
    return true;
}
