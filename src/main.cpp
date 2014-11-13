//
//  main.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/2/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

/********************** include files *****************************************/
#include <iostream>
#include <unistd.h>
#include "test.h"
#include "environment/environment.h"
#include "job/job.h"
using namespace std;
using namespace para;

void MonteCarlo(const Job &);
void Dyson(const Job &);
int main(int argc, const char *argv[])
{
    LOGGER_CONF("test.log", "test", Logger::file_on | Logger::screen_on, INFO, INFO);
    RunTest();

    string InputFile = "infile/_in_MC_2";
    para::Job Job(InputFile);
    LOGGER_CONF(ToString(Job.PID) + ".log", Job.Type, Logger::file_on | Logger::screen_on, INFO, INFO);

    if (Job.Type == "MC") {
        MonteCarlo(Job);
    }
    else if (Job.Type == "DYSON") {
        Dyson(Job);
    }
    return 0;
}

void MonteCarlo(const para::Job &Job)
{
    EnvMonteCarlo PaddyField(Job.PID);
    if (Job.DoesLoad)
        PaddyField.Load();
    else
        PaddyField.BuildNew(Job.InputFile, Job.StartFromBare);

    auto &Para = PaddyField.Para;
    auto &Grasshopper = PaddyField.Grasshopper;
    auto &Scarecrow = PaddyField.Scarecrow;

    while (Para.Counter < 10) {
        Para.Counter++;
        Grasshopper.Hop(0);

        Scarecrow.Measure();

        if (Para.Counter % 10 == 0) {
            //            Env.AddStatistics();
            Scarecrow.ReWeightEachOrder();
        }
        if (Para.Counter % 100000 == 0) {
            PaddyField.Save();
        }
        if (Para.Counter % 20 == 0) {
            PaddyField.ListenMessage();
        }
    }
}

void Dyson(const para::Job &Job)
{
    EnvDyson env(Job.PID);
    if (Job.DoesLoad) {
        env.Load();
        env.UpdateWeight();
    }
    else
        env.BuildNew(Job.InputFile, Job.StartFromBare);
    auto &Para = env.Para;

    do {
        LOG_INFO("Start Dyson Version " << Para.Version << "...")
        //TODO: do dyson

        Para.Version++;
        env.Save();
        LOG_INFO("Version " << Para.Version << " is done!")

        if (Job.DoesLoad) {
            LOG_INFO("Go to Sleep for " << Para.SleepTime << " sec...");
            sleep(Para.SleepTime);
            env.Load(); //Load new Sigma, Polar, and G, W, para will be loaded as well
            env.UpdateWeight();
        }
    } while (Job.DoesLoad);
}
