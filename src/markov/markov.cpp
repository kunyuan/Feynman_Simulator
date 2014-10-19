//
//  markov.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/17/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "markov.h"

Markov::Markov(EnvMoneCarlo *Env)
{
    _Env=Env;
    _Diag=&_Env->Diag;
    _RNG=&_Env->RNG;
}

/**
*  \brief let the Grasshopper hops for Steps
*
*  @param Steps
*/
void Markov::Hop(int && Steps)
{
    const double W1=1.0;
    const double W2=1.0;
    const double W=W1+W2;
    double x=_RNG->urn();
    if(x<W1/W)
        CreateWorm();
    else if(x<(W1+W2)/W)
        DeleteWorm();
}

void Markov::CreateWorm()
{
    
}

void Markov::DeleteWorm()
{
    
}