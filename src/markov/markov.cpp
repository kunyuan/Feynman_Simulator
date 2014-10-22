//
//  markov.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/17/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "markov.h"
#include "weight.h"

Markov::Markov(EnvMonteCarlo *Env)
{
    Beta=Env->Beta;
    Lat=&Env->Lat;
    OrderWeight=Env->OrderWeight;
    Diag = &Env->Diag;
    Worm = &Env->Diag.Worm;
    Sigma=&Env->Sigma;
    Polar=&Env->Polar;
    G=&Env->G;
    W=&Env->W;
    WormWeight=&Env->Worm;
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
    double x = RNG.urn();
    if (x < W1 / W)
        CreateWorm(W1, W2);
    else if (x < (W1+W2) / W)
        DeleteWorm(W2, W1);
}

void Markov::CreateWorm(real pcall1, real pcall2)
{
    if(Worm->Exist) return;
    WLine& w=Diag->RandomPickW();
    Vertex& vin=Diag->NeighVer(w, IN);
    Vertex& vout=Diag->NeighVer(w, OUT);
    
    int k = RandomPickK();
    int dspin = RandomPickdSpin();
    if(Diag->CanNotMoveWorm(dspin, vin) && Diag->CanNotMoveWorm(-dspin, vout)) return;
    
    Complex wWeight= W->Weight(Lat->Distance(vin.R, vout.R), vout.Tau-vin.Tau, vin.Spin, vout.Spin, true);
    
    Complex weightratio = wWeight/w.Weight;
    Complex sgn=phase(weightratio);
    
    real prob=mod(weightratio);
    real wormWeight = WormWeight->Weight(Lat->Distance(vin.R, vout.R), vout.Tau-vin.Tau);
    
    prob *= pcall2/pcall1*wormWeight*Diag->Order*2.0;
    
    if(RNG.urn()<prob)
    {
        Diag->Phase *= sgn;
        Diag->Weight *= weightratio;
        
        Worm->Exist = true;
        Worm->Ira = vin.Name;
        Worm->Masha = vout.Name;
        Worm->dSpin = dspin;
        Worm->K = k;
        
        w.IsWorm = true;
        w.Weight = wWeight;  //has to be after Diag->Weight *= ...
    }
}

void Markov::DeleteWorm(real pcall1, real pcall2)
{
    if(!Worm->Exist) return;
}