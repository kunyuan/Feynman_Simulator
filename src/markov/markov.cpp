//
//  markov.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/17/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "markov.h"
#include "weight.h"

const int MAX_K = 10000;

int RandomPickK();
int RandomPickDeltaSpin();
bool CanNotMoveWorm(int dspin, const Vertex &v);

Markov::Markov(EnvMonteCarlo *Env)
{
    Beta = Env->Beta;
    Lat = &Env->Lat;
    OrderWeight = Env->OrderWeight;
    Diag = &Env->Diag;
    Worm = &Env->Diag.Worm;
    Sigma = &Env->Sigma;
    Polar = &Env->Polar;
    G = &Env->G;
    W = &Env->W;
    WormWeight = &Env->Worm;
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
    else if (x < (W1 + W2) / W)
        DeleteWorm(W2, W1);
}
void Markov::CreateWorm(real pcall1, real pcall2)
{
    if (Worm->Exist)
        return;
    WLine &w = Diag->RandomPickW();
    Vertex &vin = Diag->NeighVer(w, IN);
    Vertex &vout = Diag->NeighVer(w, OUT);

    int k = RandomPickK();
    int dspin = RandomPickDeltaSpin();
    if (CanNotMoveWorm(dspin, vin) && CanNotMoveWorm(-dspin, vout))
        return;

    Complex wWeight = W->Weight(Lat->Distance(vin.R, vout.R), vout.Tau - vin.Tau, vin.Spin, vout.Spin, true);
    Complex weightRatio = wWeight / w.Weight;
    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    real wormWeight = WormWeight->Weight(Lat->Distance(vin.R, vout.R), vout.Tau - vin.Tau);

    prob *= pcall2 / pcall1 * wormWeight * Diag->Order * 2.0;

    if (prob >= 1.0 || RNG.urn() < prob) {
        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;

        Worm->Exist = true;
        Worm->Ira = vin.Name;
        Worm->Masha = vout.Name;
        Worm->dSpin = dspin;
        Worm->K = k;

        w.IsWorm = true;
        w.Weight = wWeight;
    }
}

void Markov::DeleteWorm(real pcall1, real pcall2)
{
    if (!Worm->Exist)
        return;
}

bool CanNotMoveWorm(int dspin, const Vertex &v)
{
    if (dspin == 1 && v.Spin[IN] == DOWN && v.Spin[OUT] == UP)
        return true;
    if (dspin == -1 && v.Spin[IN] == UP && v.Spin[OUT] == DOWN)
        return true;
    return false;
}

int RandomPickK()
{
    return RNG.irn(-MAX_K, MAX_K);
}

int RandomPickDeltaSpin()
{
    return RNG.irn(0, 2) * 2 - 1;
}
