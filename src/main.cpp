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
#include "markov.h"
using namespace std;

void MonteCarlo(Jobs &Job);
void Dyson(Jobs &Job);
int main(int argc, const char *argv[])
{
    //initialize LOGGER
    LOGGER_CONF("", "MC", Logger::file_on | Logger::screen_on, INFO, INFO);

    InitGlobalUtility();
    RunTest();

    Jobs Job;
    Job.Read();

    InitEveryOneNeedsIt();
    switch (Job.Type) {
        case MC:
            MonteCarlo(Job);
        case DYSON:
            Dyson(Job);
    }
    return 0;
}

void MonteCarlo(Jobs &Job)
{
    EnvMoneCarlo Env(Job);
    Markov GrassHopper(&Env);
    while (Env.Counter < 10000) {
        GrassHopper.Hop(10);

        Env.Measure();

        if (Env.Counter % 10) {
            //            Env.AddStatistics();
            Env.ReWeightEachOrder();
        }
        if (Env.Counter % 20) {
            Env.SaveState();
        }
        if (Env.Counter % 20) {
            Env.Annealing();
        }
    }
}

void Dyson(Jobs &Job)
{
    EnvDyson Env(Job);
}
