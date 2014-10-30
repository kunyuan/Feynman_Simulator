//
//  envMonteCarlo.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/18/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "environment.h"

using namespace std;
EnvMonteCarlo::EnvMonteCarlo()
{
    Sigma = nullptr;
    Polar = nullptr;
    G = nullptr;
    W = nullptr;
    WormWeight = nullptr;
}

EnvMonteCarlo::~EnvMonteCarlo()
{
    delete Sigma;
    delete Polar;
    delete G;
    delete W;
    delete WormWeight;
}

bool EnvMonteCarlo::BuildFromFile(string InputFile)
{
    Environment::BuildFromFile(InputFile);

    //Read more parameters
    GET(_para, Toss);
    GET(_para, Sample);
    GET(_para, Sweep);
    GET(_para, Seed);
    GET(_para, ReadFile);
    GET(_para, WormSpaceReweight);
    ReadOrderWeight();

    //Initialize utilies for MC simulations
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
    cEstimator.AddEstimator("2");
    rEstimator.AddEstimator("2");
    rEstimator.AddEstimator("3");
    return true;
}

void EnvMonteCarlo::SaveState()
{
    Environment::SaveState();
}

void EnvMonteCarlo::ReadOrderWeight(char sep)
{
    stringstream ss(_para.front());
    char sepchar;
    ss >> OrderWeight[0];
    for (int i = 1; i < Order; i++) {
        ss >> sepchar;
        if (sepchar != sep)
            ABORT("Sep char " << sepchar << " is not expected as the separator. I will expect " << sep);
        ss >> OrderWeight[i];
    }
    if (ss.fail())
        ABORT("Read OrderWeight fails!");
    _para.pop_front();
}
