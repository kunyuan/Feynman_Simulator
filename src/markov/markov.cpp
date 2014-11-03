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
    //TODO: shouldn't initialize G, W here
    G->InitializeState();
    W->InitializeState();
    WormWeight = &weight.WormWeight;

    for (int i = 0; i < NUpdates; i++)
        ProbofCall[i] = 1.0 / real(NUpdates);
    return true;
}

/**
*  \brief let the Grasshopper hops for Steps
*
*  @param Steps
*/
void Markov::Hop(int sweep)
{
    for (int i = 0; i < sweep; i++) {
        double x = RNG.urn();
        if (x < ProbofCall[0])
            CreateWorm();
        else if (x < ProbofCall[0] + ProbofCall[1])
            DeleteWorm();
        else if (x < ProbofCall[0] + ProbofCall[1] + ProbofCall[2])
            MoveWormOnG();
        else if (x < ProbofCall[0] + ProbofCall[1] + ProbofCall[2] + ProbofCall[3])
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
    vertex vin = Diag->NeighVer(w, IN);
    vertex vout = Diag->NeighVer(w, OUT);

    int k = RandomPickK();
    int dspin = RandomPickDeltaSpin();
    if (CanNotMoveWorm(dspin, Diag->Spin(vin, IN), Diag->Spin(vin, OUT)) && CanNotMoveWorm(-dspin, Diag->Spin(vout, IN), Diag->Spin(vout, OUT)))
        return;

    Complex wWeight;
    if (w->IsDelta)
        wWeight = W->WeightOfDelta(vin->R, vout->R, vin->Spin, vout->Spin, true);
    else
        wWeight = W->Weight(vin->R, vout->R, vin->Tau, vout->Tau, vin->Spin, vout->Spin, true, w->IsMeasure);
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
    wLine w = Diag->NeighW(Worm->Ira);
    if (w != Diag->NeighW(Worm->Masha))
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

    Complex wWeight;
    if (w->IsDelta)
        wWeight = W->WeightOfDelta(vin->R, vout->R, vin->Spin, vout->Spin, false);
    else
        wWeight = W->Weight(vin->R, vout->R, vin->Tau, vout->Tau, vin->Spin, vout->Spin, false, w->IsMeasure);

    Complex weightRatio = wWeight / w->Weight;
    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    prob *= ProbofCall[0] / (ProbofCall[1] * Worm->Weight * Diag->Order * 2.0);

    if (prob >= 1.0 || RNG.urn() < prob) {
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
    int dir = RandomPickDir();
    vertex vi = Worm->Ira;
    gLine g = Diag->NeighG(vi, dir);
    vertex vj = Diag->NeighVer(g, dir);
    if (vj == Worm->Masha)
        return;
    if (CanNotMoveWorm(Worm->dSpin, Diag->Spin(g), dir))
        return;

    wLine wi = Diag->NeighW(vi);
    vertex vWi = Diag->NeighVer(wi, FlipDir(vi->Dir));
    bool isWormWi;
    if (vWi == vj)
        isWormWi = true;
    else
        isWormWi = Diag->IsWorm(vWi);

    spin spinVi[2] = {vi->Spin[0], vi->Spin[1]};
    spinVi[dir] = FlipSpin(spinVi[dir]);

    Complex wiWeight;
    if (wi->IsDelta)
        wiWeight = W->WeightOfDelta(vi->Dir, vi->R, vWi->R, spinVi, vWi->Spin, isWormWi);
    else
        wiWeight = W->Weight(vi->Dir, vi->R, vWi->R, vi->Tau, vWi->Tau, spinVi, vWi->Spin, isWormWi, wi->IsMeasure);

    wLine wj = Diag->NeighW(vj);
    vertex vWj = Diag->NeighVer(wj, FlipDir(vj->Dir));

    spin spinVj[2] = {vj->Spin[0], vj->Spin[1]};
    spinVj[FlipDir(dir)] = FlipSpin(spinVj[FlipDir(dir)]);

    Complex wjWeight;
    if (wj->IsDelta)
        wjWeight = W->WeightOfDelta(vj->Dir, vj->R, vWj->R, spinVj, vWj->Spin, true);
    else
        wjWeight = W->Weight(vj->Dir, vj->R, vWj->R, vj->Tau, vWj->Tau, spinVj, vWj->Spin, true, wj->IsMeasure);
}

void Markov::MoveWormOnW()
{
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

int RandomPickK()
{
    return RNG.irn(-MAX_K, MAX_K - 1);
}

int RandomPickDeltaSpin()
{
    return RNG.irn(0, 1) * 2 - 1;
}

int RandomPickDir()
{
    return RNG.irn(0, 1);
}
