//
//  markov.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/17/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "markov.h"

Markov::Markov(EnvMonteCarlo *Env)
{
    Beta=Env->Beta;
    Lat=&Env->Lat;
    OrderWeight=Env->OrderWeight;
    Diag = &Env->Diag;
    RNG = &Env->RNG;
    Sigma=&Env->Sigma;
}

/**
*  \brief let the Grasshopper hops for Steps
*
*  @param Steps
*/
void Markov::Hop(int &&Steps)
{
    const double W1 = 1.0;
    const double W2 = 1.0;
    const double W = W1 + W2;
    double x = RNG->urn();
    if (x < W1 / W)
        CreateWorm();
    else if (x < (W1 + W2) / W)
        DeleteWorm();
}

void Markov::CreateWorm()
{
}

void Markov::DeleteWorm()
{
}