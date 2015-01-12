//
//  envMonteCarlo.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/18/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "environment.h"
#include "utility/dictionary.h"

using namespace std;
using namespace para;

const string ParaKey = "Para";
const string ConfigKey = "Config";
const string WeightKey = "Weight";
const string HistKey = "Histogram";
const string EstimatorsKey = "Estimators";

EnvMonteCarlo::EnvMonteCarlo(const para::Job &job, bool IsAllTauSymmetric)
    : Job(job), Weight(IsAllTauSymmetric)
{
}

bool EnvMonteCarlo::BuildNew()
{
    LOGGER_CONF(Job.LogFile, Job.Type, Logger::file_on | Logger::screen_on, INFO, INFO);
    //Read more stuff for the state of MC only
    Dictionary para_;
    para_.Load(Job.InputFile);
    Para.FromDict(para_.Get<Dictionary>(ParaKey));
    //Load GW weight from a global file shared by other MC processes
    Dictionary GW_;
    GW_.BigLoad(Job.WeightFile);
    Weight.FromDict(GW_, weight::GW, Para);
    
    //TEST//
    Weight.SetDiagCounter(Para);
    
    Weight.BuildNew(weight::SigmaPolar, Para);
    Diag.BuildNew(Para.Lat, *Weight.G, *Weight.W);
    Grasshopper.BuildNew(Para, Diag, Weight);
    Scarecrow.BuildNew(Para, Diag, Weight);
    para_[ConfigKey] = Diag.ToDict();
    para_.Save(Job.ParaFile, "w");
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
    Dictionary para_;
    para_.Load(Job.ParaFile);
    Para.FromDict(para_.Get<Dictionary>(ParaKey));
    Dictionary statis_;
    statis_.BigLoad(Job.StatisticsFile);
    Weight.FromDict(statis_, weight::GW, Para);
    Weight.FromDict(statis_, weight::SigmaPolar, Para);
    Diag.FromDict(para_.Get<Dictionary>(ConfigKey), Para.Lat, *Weight.G, *Weight.W);
    Scarecrow.FromDict(statis_, Para, Diag, Weight);
    Grasshopper.BuildNew(Para, Diag, Weight);
    return true;
}

void EnvMonteCarlo::Save()
{
    Dictionary para_;
    para_[ParaKey] = Para.ToDict();
    para_[ConfigKey] = Diag.ToDict();
    para_["PID"] = Job.PID;
    para_.Save(Job.ParaFile, "w");
    Dictionary statis_ = Weight.ToDict(weight::GW | weight::SigmaPolar);
    statis_.Update(Scarecrow.ToDict());
    statis_.BigSave(Job.StatisticsFile);
}
void EnvMonteCarlo::DeleteSavedFiles()
{
    system(("rm " + Job.ParaFile).c_str());
    system(("rm " + Job.StatisticsFile).c_str());
    system(("rm " + Job.WeightFile).c_str());
}

/**
*  Adjust everything according to new parameters, like new Beta, Jcp
*/
bool EnvMonteCarlo::ListenToMessage()
{
    LOG_INFO("Start reweighting...");
    Message Message_;
    if (!Message_.Load(Job.MessageFile))
        return false;
    if (Para.Version >= Message_.Version) {
        LOG_INFO("Status has not been updated yet since the last reweighting!");
        return false;
    }
    Para.UpdateWithMessage(Message_);
    Dictionary weight_;
    weight_.Load(Job.WeightFile);
    Weight.FromDict(weight_, weight::GW, Para);
    Weight.ReWeight(weight::GW | weight::SigmaPolar, Para);
    Grasshopper.ReWeight(Para);
    Scarecrow.ReWeight();
    LOG_INFO("Reweighted to:\n" << Message_.PrettyString());
    return true;
}
