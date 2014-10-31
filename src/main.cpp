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
#include "markov/markov.h"
#include "markov/markov_monitor.h"
using namespace std;

void MonteCarlo(EnvMonteCarlo &);
void Dyson(EnvDyson &);
int main(int argc, const char *argv[])
{
    //initialize LOGGER
    LOGGER_CONF("test.log", "TEST", Logger::file_on | Logger::screen_on, INFO, INFO);

    InitGlobalUtility();
    RunTest();
    //    InitEveryOneNeedsIt();
    string InputFile = "infile/_in_MC_1";

    switch (GetJobsType(InputFile)) {
        case MC: {
            EnvMonteCarlo env;
            env.BuildFromFile(InputFile);
            MonteCarlo(env);
        }
        case DYSON:
            EnvDyson env;
            env.BuildFromFile(InputFile);
            Dyson(env);
    }
    return 0;
}

void MonteCarlo(EnvMonteCarlo &PaddyField)
{
    Markov GrassHopper(&PaddyField);
    MarkovMonitor Scarecrow(&PaddyField);
    while (PaddyField.Counter < 10) {
        PaddyField.Counter++;
        //        GrassHopper.Hop(10);

        Scarecrow.Measure();

        if (PaddyField.Counter % 10) {
            //            Env.AddStatistics();
            Scarecrow.ReWeightEachOrder();
        }
        if (PaddyField.Counter % 100000) {
            PaddyField.SaveState();
        }
        if (PaddyField.Counter % 20) {
            Scarecrow.Annealing();
        }
    }
}

void Dyson(EnvDyson &env)
{
}
