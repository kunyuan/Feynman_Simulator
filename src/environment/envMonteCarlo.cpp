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

EnvMonteCarlo::EnvMonteCarlo(const para::Job& job, bool IsAllTauSymmetric)
    : Job(job)
    , Weight(IsAllTauSymmetric)
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
    if (Para.runGamma3)
        Weight.FromDict(GW_, weight::GammaGW, Para);

    //    Weight.SetDiagCounter(Para);//Test for DiagCounter
    //    Weight.SetTest(Para);//Test for WeightTest

    Weight.BuildNew(weight::SigmaPolar, Para);
    Diag.BuildNew(Para.Lat, *Weight.G, *Weight.W);
    Markov.BuildNew(Para, Diag, Weight);
    MarkovMonitor.BuildNew(Para, Diag, Weight);
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
    LOG_INFO("Trying to load " << Job.ParaFile);
    bool DoesParaFileExit;
    try {
        para_.Load(Job.ParaFile);
        DoesParaFileExit = true;
    }
    catch (IOInvalid e) {
        LOG_WARNING("Load " << Job.ParaFile << " failed, use " << Job.InputFile << " instead!");
        para_.Load(Job.InputFile);
        DoesParaFileExit = false;
    }
    Para.FromDict(para_.Get<Dictionary>(ParaKey));
    Dictionary GW_;
    GW_.BigLoad(Job.WeightFile);
    Dictionary statis_;
    statis_.BigLoad(Job.StatisticsFile);
    Weight.FromDict(GW_, weight::GW, Para);
    Weight.FromDict(statis_, weight::SigmaPolar, Para);
    if(Para.runGamma3)
    {
        Weight.FromDict(GW_, weight::GammaGW, Para);
        Weight.FromDict(statis_, weight::GammaGWStatis, Para);
    }
    LOG_INFO(DoesParaFileExit);
    if (DoesParaFileExit)
        Diag.FromDict(para_.Get<Dictionary>(ConfigKey), Para.Lat, *Weight.G, *Weight.W);
    else
        Diag.BuildNew(Para.Lat, *Weight.G, *Weight.W);
    MarkovMonitor.FromDict(statis_, Para, Diag, Weight);
    Markov.BuildNew(Para, Diag, Weight);
    return true;
}

void EnvMonteCarlo::Save()
{
    if (Diag.HasGammaGW==0) {
        LOG_INFO("Start saving data...");
        Dictionary para_;
        para_[ParaKey] = Para.ToDict();
        para_[ConfigKey] = Diag.ToDict();
        para_["PID"] = Job.PID;
        para_.Save(Job.ParaFile, "w");
        Dictionary statis_;
        if (Para.runGamma3)
            statis_ = Weight.ToDict(weight::SigmaPolar | weight::GammaGWStatis);
        else
            statis_ = Weight.ToDict(weight::SigmaPolar);
        statis_.Update(MarkovMonitor.ToDict());
        statis_.BigSave(Job.StatisticsFile);
        LOG_INFO("Saving data is done!");
    }
}

void EnvMonteCarlo::DeleteSavedFiles()
{
    system(("rm " + Job.ParaFile).c_str());
    system(("rm " + Job.StatisticsFile).c_str());
    system(("rm " + Job.WeightFile).c_str());
}

void EnvMonteCarlo::AdjustOrderReWeight()
{
    LOG_INFO("Start adjusting OrderReweight...");
    if (MarkovMonitor.AdjustOrderReWeight()) {
        Markov.Reset(Para, Diag, Weight);
        string str;
        for (int i = 0; i <= Para.Order; i++)
            str += ToString((Para.OrderReWeight[i])) + "  ";
        if(Para.runGamma3){
            LOG_INFO("Reweighted to:\n" + str + "\nWorm Reweighted to:\n" + ToString(Para.WormSpaceReweight)
                     + "\nPolar Reweighted to:\n" + ToString(Para.PolarReweight)
                     + "\nGammaG Reweighted to:\n" + ToString(Para.GammaGReweight)
                     + "\nGammaW Reweighted to:\n" + ToString(Para.GammaWReweight));
        }else{
            LOG_INFO("Reweighted to:\n" + str + "\nWorm Reweighted to:\n" + ToString(Para.WormSpaceReweight)
                     + "\nPolar Reweighted to:\n" + ToString(Para.PolarReweight));
        }
    }
    else {
        string str;
        for (int i = 0; i <= Para.Order; i++)
            str += ToString((MarkovMonitor.PhyEstimator[i].Norm())) + "  ";
        LOG_INFO("Number of samples is too small, adjust later.\n"
                 << "Norm of different orders: " << str);
    }
}
/**
*  Adjust everything according to new parameters, like new Beta, Jcp
*/
bool EnvMonteCarlo::ListenToMessage()
{
    if(!Para.runGamma3 || Diag.HasGammaGW==0){
        if (Para.runGamma3){
            LOG_WARNING("Beginning:");
            LOG_WARNING(Diag.HasGammaGW);
        }

        LOG_INFO("Start Annealing...");
        Message Message_;
        if (!Message_.Load(Job.MessageFile))
            return false;
        if (Para.Version >= Message_.Version) {
            LOG_INFO("Status has not been updated yet since the last annealing!");
            return false;
        }
        Dictionary weight_;
        try {
            weight_.BigLoad(Job.WeightFile);
        }
        catch (IOInvalid e) {
            LOG_WARNING("Annealing Failed!");
            return false;
        }
        Para.UpdateWithMessage(Message_);

        Weight.FromDict(weight_, weight::GW, Para);

        if (Para.runGamma3){
            Weight.FromDict(weight_, weight::GammaGW, Para);
        }

        Weight.Anneal(Para);
        LOG_WARNING("Calling Diagram reset!");

        if (Para.runGamma3) {
            LOG_WARNING("After:");
            LOG_WARNING(Diag.HasGammaGW);
        }

        Diag.Reset(Para.Lat, *Weight.G, *Weight.W);
        Markov.Reset(Para, Diag, Weight);
        MarkovMonitor.Reset(Para, Diag, Weight);
        MarkovMonitor.SqueezeStatistics(Message_.SqueezeFactor);
        LOG_INFO("Annealled to " << Message_.PrettyString()
                                 << "\nwith squeeze factor" << Message_.SqueezeFactor);
        if(Diag.CheckDiagram()){
            LOG_INFO("Diagram Check pass.");
        }else {
            ABORT("Diagram Check didn't pass when annealing!");
        }

        return true;
    }
}
