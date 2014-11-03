//
//  markov.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/17/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "markov.h"

const int MAX_K = 10000;

int RandomPickK();
int RandomPickDeltaSpin();
bool CanNotMoveWorm(int dspin, vertex v);

bool Markov::BuildNew(ParameterMC &para, Diagram &diag, weight::Weight &weight)
{
    Beta = para.Beta;
    Lat = &para.Lat;
    OrderWeight = para.OrderReWeight;
    Diag = &diag;
    Worm = &diag.Worm;
    Sigma = weight.Sigma;
    Polar = weight.Polar;
    G = weight.G;
    W = weight.W;
    //TODO: shouldn't initialize G, W here
    G->InitializeState();
    W->InitializeState();
    WormWeight = weight.WormWeight;

    ProbofCall[0] = 0.50;
    ProbofCall[1] = 0.50;
    return true;
}

/**
*  \brief let the Grasshopper hops for Steps
*
*  @param Steps
*/
void Markov::Hop(int sweep)
{
    double x = RNG.urn();
    if (x < ProbofCall[0])
        CreateWorm();
    else if (x < ProbofCall[0] + ProbofCall[1])
        DeleteWorm();
}

void Markov::CreateWorm()
{
    if (Worm->Exist)
        return;

    wLine w = Diag->RandomPickW();
    vertex vin = Diag->NeighVer(w, IN);
    vertex vout = Diag->NeighVer(w, OUT);

    int k = RandomPickK();
    int dspin = RandomPickDeltaSpin();
    if (CanNotMoveWorm(dspin, vin) && CanNotMoveWorm(-dspin, vout))
        return;

    Complex wWeight = W->Weight(vin->R, vout->R, vin->Tau, vout->Tau, vin->Spin, vout->Spin, true);
    Complex weightRatio = wWeight / w->Weight;
    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    real wormWeight = WormWeight->Weight(vin->R, vout->R, vin->Tau, vout->Tau);

    prob *= ProbofCall[1] / ProbofCall[0] * wormWeight * Diag->Order * 2.0;

    if (prob >= 1.0 || RNG.urn() < prob) {
        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;

        Worm->Exist = true;
        Worm->Ira = vin;
        Worm->Masha = vout;
        Worm->dSpin = dspin;
        Worm->K = k;

        w->IsWorm = true;
        w->K -= k;
        w->Weight = wWeight;
    }
}

void Markov::DeleteWorm()
{
    if (!Worm->Exist)
        return;
    if (Diag->NeighW(Worm->Ira) != Diag->NeighW(Worm->Masha))
        return;
}

bool CanNotMoveWorm(int dspin, vertex v)
{
    if (dspin == 1 && v->Spin[IN] == DOWN && v->Spin[OUT] == UP)
        return true;
    if (dspin == -1 && v->Spin[IN] == UP && v->Spin[OUT] == DOWN)
        return true;
    return false;
}

int RandomPickK()
{
    return RNG.irn(-MAX_K, MAX_K - 1);
}

int RandomPickDeltaSpin()
{
    return RNG.irn(0, 1) * 2 - 1;
}
