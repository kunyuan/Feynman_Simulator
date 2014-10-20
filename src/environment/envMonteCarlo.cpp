//
//  envMonteCarlo.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/18/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "environment.h"

EnvMoneCarlo::EnvMoneCarlo(Jobs &Job)
    : Environment(Job),
      //The base class is constructed first!!!
      Sigma(Lat, Beta, StateFile),
      Pi(Lat, Beta, StateFile)
{
    //Initialize random number generator
    RNG.Reset(Job.Seed);

    //Add estimators
    cEstimator.AddEstimator("1");
    rEstimator.AddEstimator("1");
}

void EnvMoneCarlo::Annealing()
{
    SqueezeStatistics();
}

void EnvMoneCarlo::SqueezeStatistics()
{
}

void EnvMoneCarlo::ReWeightEachOrder()
{
}

void EnvMoneCarlo::Measure()
{
    //    cEstimator[0].Measure(<#const Complex &#>)
}

void EnvMoneCarlo::AddStatistics()
{
    cEstimator.AddStatistics();
    rEstimator.AddStatistics();
}

void EnvMoneCarlo::SaveState()
{
    Environment::SaveState();
}
