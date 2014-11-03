//
//  main.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/2/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

/********************** include files *****************************************/
#include <iostream>
#include "definition.h"
#include "initialization.h"
#include "test.h"
#include "environment/environment.h"
#include "parameter/job.h"
#include "markov/markov.h"
#include "markov/markov_monitor.h"
using namespace std;

void MonteCarlo(EnvMonteCarlo &);
void Dyson(EnvDyson &);
int main(int argc, const char *argv[])
{
    //initialize LOGGER
    LOGGER_CONF("test.log", "TEST", Logger::file_on | Logger::screen_on, INFO, INFO);

    RunTest();
    //    InitEveryOneNeedsIt();
    string InputFile = "infile/_in_MC_1";
    job Job(InputFile);
    LOGGER_CONF(ToString(Job.PID) + ".log", Job.Type, Logger::file_on | Logger::screen_on, INFO, INFO);

    if (Job.Type == "MC") {
        EnvMonteCarlo env(Job.PID);
        env.BuildNew(InputFile);
        MonteCarlo(env);
    }
    else if (Job.Type == "DYSON") {
        EnvDyson env(Job.PID);
        env.BuildNew(InputFile);
        Dyson(env);
    }
    return 0;
}

void MonteCarlo(EnvMonteCarlo &PaddyField)
{

    while (PaddyField.Para.Counter < 10) {
        PaddyField.Para.Counter++;
        //        GrassHopper.Hop(10);

        PaddyField.Scarecrow.Measure();

        if (PaddyField.Para.Counter % 10) {
            //            Env.AddStatistics();
            PaddyField.Scarecrow.ReWeightEachOrder();
        }
        if (PaddyField.Para.Counter % 100000) {
            PaddyField.Save();
        }
        if (PaddyField.Para.Counter % 20) {
            PaddyField.Scarecrow.Annealing();
        }
    }
}

void Dyson(EnvDyson &env)
{
}
