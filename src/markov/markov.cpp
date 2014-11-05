//
//  markov.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/17/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "markov.h"
#include "math.h"

const int MAX_K = 10000;

int RandomPickK();
int RandomPickDir();
int RandomPickDeltaSpin();
bool CanNotMoveWorm(int dspin, spin sin, spin sout);
bool CanNotMoveWorm(int dspin, spin sin, int dir);

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
    WormWeight = &weight.WormWeight;
    RNG = &para.RNG;

    for (int i = 0; i < NUpdates; i++) {
        ProbofCall[i] = 1.0 / real(NUpdates);
        for (int j = i; j < NUpdates; j++)
            SumofProbofCall[j] += ProbofCall[i];
    }
    return true;
}

void Markov::ReWeight(ParameterMC &para)
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
    vertex vi = Worm->Ira;
    int dir = RandomPickDir();
    gLine g = vi->NeighG(dir);
    vertex vj = g->NeighVer(dir);
    if (vj == Worm->Masha || CanNotMoveWorm(Worm->dSpin, g->Spin(), dir))
        return;

    wLine wi = vi->NeighW();
    vertex vWi = wi->NeighVer(INVERSE(vi->Dir));
    bool isWormWi;
    if (vWi == vj)
        isWormWi = true;
    else
        isWormWi = Diag->IsWorm(vWi);

    spin spinVi[2] = {vi->Spin(0), vi->Spin(1)};
    spinVi[dir] = FLIP(spinVi[dir]);

    Complex wiWeight = W->Weight(vi->Dir, vi->R, vWi->R, vi->Tau, vWi->Tau,
                             spinVi, vWi->Spin(), isWormWi, wi->IsMeasure, wi->IsDelta);

    wLine wj = vj->NeighW();
    vertex vWj = wj->NeighVer(INVERSE(vj->Dir));

    spin spinVj[2] = {vj->Spin(0), vj->Spin(1)};
    spinVj[INVERSE(dir)] = FLIP(spinVj[INVERSE(dir)]);

    Complex wjWeight = W->Weight(vj->Dir, vj->R, vWj->R, vj->Tau, vWj->Tau,
                                 spinVj, vWj->Spin(), true, wj->IsMeasure, wi->IsDelta);

    Complex gWeight = G->Weight(INVERSE(dir), vi->R, vj->R, vi->Tau, vj->Tau,
                                spinVi[dir], spinVj[INVERSE(dir)], g->IsMeasure);

    Complex weightRatio = wiWeight * wjWeight * gWeight / (g->Weight * wi->Weight * wj->Weight);
    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    real wormWeight = WormWeight->Weight(vj->R, Worm->Masha->R, vj->Tau, Worm->Masha->Tau);

    prob *= wormWeight / Worm->Weight;

    if (prob >= 1.0 || RNG->urn() < prob) {
        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;

        g->Weight = gWeight;
        g->K = g->K - pow(-1, dir) * Worm->K;

        wi->Weight = wiWeight;
        wi->IsWorm = isWormWi;

        wj->Weight = wjWeight;
        wj->IsWorm = true;

        vi->SetSpin(spinVi);
        vj->SetSpin(spinVj);

        Worm->Ira = vj;
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
