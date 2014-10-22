//
//  envMonteCarlo.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/18/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "environment.h"

EnvMonteCarlo::EnvMonteCarlo(Jobs &Job)
    : Environment(Job),
      //The base class is constructed first!!!
      Sigma(Lat, Beta, Order),
      Polar(Lat, Beta, Order),
      G(Lat, Beta, Order),
      W(Lat, Beta, Order)
{
    //Initialize random number generator
    RNG.Reset(Job.Seed);

    //Add estimators
    cEstimator.AddEstimator("1");
    rEstimator.AddEstimator("1");
}

void EnvMonteCarlo::SaveState()
{
    Environment::SaveState();
}
