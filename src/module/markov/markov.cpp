//
//  markov.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/17/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "markov.h"
#include "math.h"
#include "module/weight/weight.h"
#include "module/diagram/diagram.h"
#include "module/parameter/parameter.h"

using namespace std;
using namespace diag;
using namespace para;
using namespace mc;

const int MAX_K = 10000;

int RandomPickK();
int RandomPickDir();
int RandomPickDeltaSpin();
bool CanNotMoveWorm(int dspin, spin sin, spin sout);
bool CanNotMoveWorm(int dspin, spin sin, int dir);

bool Markov::BuildNew(ParaMC &para, Diagram &diag, weight::Weight &weight)
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
    WormWeight = &weight.WormWeight;
    RNG = &para.RNG;

    for (int i = 0; i < NUpdates; i++) {
        ProbofCall[i] = 1.0 / real(NUpdates);
        for (int j = i; j < NUpdates; j++)
            SumofProbofCall[j] += ProbofCall[i];
    }
    return true;
}

void Markov::ReWeight(ParaMC &para)
{
    Beta = para.Beta;
    OrderWeight = para.OrderReWeight;

    for (int i = 0; i < NUpdates; i++)
        ProbofCall[i] = 1.0 / real(NUpdates);
    //TODO: consider what have to be reweighted here
}

/**
*  \brief let the Grasshopper hops for Steps
*
*  @param Steps
*/
void Markov::Hop(int sweep)
{
    for (int i = 0; i < sweep; i++) {
        double x = RNG->urn();
        if (x < SumofProbofCall[0])
            CreateWorm();
        else if (x < SumofProbofCall[1])
            DeleteWorm();
        else if (x < SumofProbofCall[2])
            MoveWormOnG();
        else if (x < SumofProbofCall[3])
            MoveWormOnW();
    }
}

/**
 *  Create Ira and Masha on a wline
 */
void Markov::CreateWorm()
{
    if (Worm->Exist)
        return;

    wLine w = Diag->RandomPickW();
    vertex vin = w->NeighVer(IN);
    vertex vout = w->NeighVer(OUT);

    int k = RandomPickK();
    int dspin = RandomPickDeltaSpin();
    if (CanNotMoveWorm(dspin, vin->Spin(IN), vin->Spin(OUT)) && CanNotMoveWorm(-dspin, vout->Spin(IN), vout->Spin(OUT)))
        return;

    Complex wWeight = W->Weight(vin->R, vout->R, vin->Tau, vout->Tau,
                                vin->Spin(), vout->Spin(), true, w->IsMeasure, w->IsDelta);
    Complex weightRatio = wWeight / w->Weight;
    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    real wormWeight = WormWeight->Weight(vin->R, vout->R, vin->Tau, vout->Tau);

    prob *= ProbofCall[1] / ProbofCall[0] * wormWeight * Diag->Order * 2.0;

    if (prob >= 1.0 || RNG->urn() < prob) {
        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;

        Worm->Exist = true;
        Worm->Ira = vin;
        Worm->Masha = vout;
        Worm->dSpin = dspin;
        Worm->K = k;
        Worm->Weight = wormWeight;

        w->IsWorm = true;
        w->K -= k;
        w->Weight = wWeight;
    }
}

/**
 *  Delete Ira and Masha on the same wline
 */
void Markov::DeleteWorm()
{
    if (!Worm->Exist)
        return;
    wLine w = Worm->Ira->NeighW();
    if (w != Worm->Masha->NeighW())
        return;
    vertex vin, vout;
    int k;
    if (Worm->Ira->Dir == IN) {
        vin = Worm->Ira;
        vout = Worm->Masha;
        k = w->K + Worm->K;
    }
    else {
        vin = Worm->Masha;
        vout = Worm->Ira;
        k = w->K - Worm->K;
    }

    Complex wWeight = W->Weight(vin->R, vout->R, vin->Tau, vout->Tau,
                                vin->Spin(), vout->Spin(), false, w->IsMeasure, w->IsDelta);

    Complex weightRatio = wWeight / w->Weight;
    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    prob *= ProbofCall[0] / (ProbofCall[1] * Worm->Weight * Diag->Order * 2.0);

    if (prob >= 1.0 || RNG->urn() < prob) {
        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;

        Worm->Exist = false;

        w->IsWorm = false;
        w->K = k;
        w->Weight = wWeight;
    }
}

/**
 *  Move Ira along a GLine
 */
void Markov::MoveWormOnG()
{
    if (!Worm->Exist)
        return;
    vertex v1 = Worm->Ira;
    int dir = RandomPickDir();
    gLine g = v1->NeighG(dir);
    vertex v2 = g->NeighVer(dir);
    if (v2 == Worm->Masha || CanNotMoveWorm(Worm->dSpin, g->Spin(), dir))
        return;

    wLine w1 = v1->NeighW();
    vertex vW1 = w1->NeighVer(INVERSE(v1->Dir));
    bool isWormW1;
    if (vW1 == v2)
        isWormW1 = true;
    else
        isWormW1 = Diag->IsWorm(vW1);

    spin spinV1[2] = {v1->Spin(0), v1->Spin(1)};
    spinV1[dir] = FLIP(spinV1[dir]);

    Complex w1Weight = W->Weight(v1->Dir, v1->R, vW1->R, v1->Tau, vW1->Tau,
                             spinV1, vW1->Spin(), isWormW1, w1->IsMeasure, w1->IsDelta);

    wLine w2 = v2->NeighW();
    vertex vW2 = w2->NeighVer(INVERSE(v2->Dir));

    spin spinV2[2] = {v2->Spin(0), v2->Spin(1)};
    spinV2[INVERSE(dir)] = FLIP(spinV2[INVERSE(dir)]);

    Complex w2Weight = W->Weight(v2->Dir, v2->R, vW2->R, v2->Tau, vW2->Tau,
                                 spinV2, vW2->Spin(), true, w2->IsMeasure, w2->IsDelta);

    Complex gWeight = G->Weight(INVERSE(dir), v1->R, v2->R, v1->Tau, v2->Tau,
                                spinV1[dir], spinV2[INVERSE(dir)], g->IsMeasure);

    Complex weightRatio = w1Weight * w2Weight * gWeight / (g->Weight * w1->Weight * w2->Weight);
    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    real wormWeight = WormWeight->Weight(v2->R, Worm->Masha->R, v2->Tau, Worm->Masha->Tau);

    prob *= wormWeight / Worm->Weight;

    if (prob >= 1.0 || RNG->urn() < prob) {
        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;

        g->Weight = gWeight;
        g->K = g->K - pow(-1, dir) * Worm->K;

        w1->Weight = w1Weight;
        w1->IsWorm = isWormW1;

        w2->Weight = w2Weight;
        w2->IsWorm = true;

        v1->SetSpin(spinV1);
        v2->SetSpin(spinV2);

        Worm->Ira = v2;
        Worm->Weight = wormWeight;
    }
}

/**
 *  Move Ira along a wline
 */
void Markov::MoveWormOnW()
{
}

int Markov::RandomPickK()
{
    return RNG->irn(-MAX_K, MAX_K - 1);
}

int Markov::RandomPickDeltaSpin()
{
    return RNG->irn(0, 1) * 2 - 1;
}

int Markov::RandomPickDir()
{
    return RNG->irn(0, 1);
}

/**
 *  determine whether a Ira can move around to another vertex
 *  used in CreateWorm
 *
 *  @param dspin the spin current of Ira
 *  @param sin   the spin incoming of vertex Ira
 *  @param sout  the spin outgoing of vertex Ira
 *
 *  @return true: cannot move; false: can move
 */
bool CanNotMoveWorm(int dspin, spin sin, spin sout)
{
    if (dspin == 1 && sin == DOWN && sout == UP)
        return true;
    if (dspin == -1 && sin == UP && sout == DOWN)
        return true;
    return false;
}

/**
 *  determine whether a Ira can move on a gline
 *  used in MoveWormOnG
 *
 *  @param dspin the spin current of Ira
 *  @param sg    the spin on Gline
 *  @param dir   the Gline is incoming or outgoing of Ira
 *
 *  @return true: cannot move; false: can move
 */
bool CanNotMoveWorm(int dspin, spin sg, int dir)
{
    if (dspin == 1) {
        if (dir == IN && sg == DOWN)
            return true;
        if (dir == OUT && sg == UP)
            return true;
    }
    else {
        if (dir == IN && sg == UP)
            return true;
        if (dir == OUT && sg == DOWN)
            return true;
    }
    return false;
}
