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
#include "parameter/job.h"
#include "markov/markov.h"
#include "markov/markov_monitor.h"
using namespace std;

void MonteCarlo(const job &);
void Dyson(const job &);
int main(int argc, const char *argv[])
{
    LOGGER_CONF("test.log", "test", Logger::file_on | Logger::screen_on, INFO, INFO);
    RunTest();

    string InputFile = "infile/_in_DYSON_1";
    job Job(InputFile);
    LOGGER_CONF(ToString(Job.PID) + ".log", Job.Type, Logger::file_on | Logger::screen_on, INFO, INFO);

    if (Job.Type == "MC") {
        MonteCarlo(Job);
    }
    else if (Job.Type == "DYSON") {
        Dyson(Job);
    }
    return 0;
}

void MonteCarlo(const job &Job)
{
    EnvMonteCarlo PaddyField(Job.PID);
    if (Job.DoesLoad)
        PaddyField.BuildNew(Job.InputFile, Job.StartFromBare);
    else
        PaddyField.Load();

    auto &Para = PaddyField.Para;
    auto &Grasshopper = PaddyField.Grasshopper;
    auto &Scarecrow = PaddyField.Scarecrow;

    while (Para.Counter < 10) {
        Para.Counter++;
        Grasshopper.Hop(0);

        Scarecrow.Measure();

        if (Para.Counter % 10) {
            //            Env.AddStatistics();
            Scarecrow.ReWeightEachOrder();
        }
        if (Para.Counter % 100000) {
            PaddyField.Save();
        }
        if (Para.Counter % 20) {
            Scarecrow.Annealing();
        }
    }
}

void Dyson(const job &Job)
{
    EnvDyson env(Job.PID);
    if (Job.DoesLoad)
        env.BuildNew(Job.InputFile, Job.StartFromBare);
    else
        env.Load();
    auto &Para = env.Para;

    while (true) {
        Para.Version++;
        LOG_INFO("Sleep");
        sleep(Para.SleepTime);
    }
}
