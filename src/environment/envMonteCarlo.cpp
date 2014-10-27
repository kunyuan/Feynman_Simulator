//
//  envMonteCarlo.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/18/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "environment.h"
#include <string>

EnvMonteCarlo::EnvMonteCarlo()
{
    OrderWeight = nullptr;
    Sigma = nullptr;
    Polar = nullptr;
    G = nullptr;
    W = nullptr;
    WormWeight = nullptr;
}

EnvMonteCarlo::~EnvMonteCarlo()
{
    delete OrderWeight;
    delete Sigma;
    delete Polar;
    delete G;
    delete W;
    delete WormWeight;
}

bool EnvMonteCarlo::BuildFromFile(string InputFile)
{
    Environment::BuildFromFile(InputFile);

    GET(_para, Toss);
    GET(_para, Sample);
    GET(_para, Sweep);
    GET(_para, Seed);
    GET(_para, WormSpaceReweight);

    Counter = 0;

    Sigma = new Weight::Sigma(*Lat, Beta, Order);
    Polar = new Weight::Polar(*Lat, Beta, Order);
    G = new Weight::G(*Lat, Beta, Order);
    W = new Weight::W(*Lat, Beta, Order);
    WormWeight = new Weight::Worm();

    Diag.Build(Lat, G, W);

    //    //Initialize random number generator
    RNG.Reset(Seed);

    //Add estimators
    cEstimator.AddEstimator("1");
    rEstimator.AddEstimator("1");
    return true;
}

void EnvMonteCarlo::SaveState()
{
    Environment::SaveState();
}
