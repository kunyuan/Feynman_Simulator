//
//  markov.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/17/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//


#include "markov.h"
#include "math.h"
#include "utility/momentum.h"
#include "module/observable/weight.h"
#include "module/diagram/diagram.h"
#include "module/parameter/parameter.h"
#include "lattice/lattice.h"

using namespace std;
using namespace diag;
using namespace para;
using namespace mc;

#define SIGN(x) ((x)==IN? 1:-1)

bool CanNotMoveWorm(int dspin, spin sin, spin sout);
bool CanNotMoveWorm(int dspin, spin sin, int dir);

bool Markov::BuildNew(ParaMC &para, Diagram &diag, weight::Weight &weight)
{
    Beta = para.Beta;
    Order = para.Order;
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
        else if (x < SumofProbofCall[4])
            Reconnect();
        else if (x < SumofProbofCall[5])
            AddInteraction();
        else if (x < SumofProbofCall[6])
            DeleteInteraction();
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

    Momentum k = RandomPickK();
    //TODO: Hash check for k
    
    int dspin = RandomPickDeltaSpin();
    if (CanNotMoveWorm(dspin, vin->Spin(IN), vin->Spin(OUT))
        && CanNotMoveWorm(-dspin, vout->Spin(IN), vout->Spin(OUT)))
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
    vertex Ira = Worm->Ira;
    vertex Masha = Worm->Masha;
    
    wLine w = Ira->NeighW();
    if (!(w == Masha->NeighW()))
        return;
    Momentum k = SIGN(Ira->Dir)* w->K + Worm->K;
    //TODO: Hash check for k

    Complex wWeight = W->Weight(Ira->Dir, Ira->R, Masha->R, Ira->Tau, Masha->Tau,
                                Ira->Spin(), Masha->Spin(), false, w->IsMeasure, w->IsDelta);

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
    
    vertex Ira = Worm->Ira;
    vertex Masha = Worm->Masha;
    
    vertex v1 = Ira;
    int dir = RandomPickDir();
    gLine g = v1->NeighG(dir);
    vertex v2 = g->NeighVer(dir);

    if (v2 == Ira || v2 == Masha || CanNotMoveWorm(Worm->dSpin, g->Spin(), dir))
        return;
    Momentum k = g->K - SIGN(dir) * Worm->K;
    //TODO: Hash check for G(knew)
    //    if(!IsKValidG(k)) return;

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

    real wormWeight = WormWeight->Weight(v2->R, Masha->R, v2->Tau, Masha->Tau);

    prob *= wormWeight / Worm->Weight;

    if (prob >= 1.0 || RNG->urn() < prob) {
        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;

        g->Weight = gWeight;
        g->K = k;

        w1->Weight = w1Weight;
        w1->IsWorm = isWormW1;

        w2->Weight = w2Weight;
        w2->IsWorm = true;

        v1->SetSpin(spinV1);
        v2->SetSpin(spinV2);

        Ira = v2;
        Worm->Weight = wormWeight;
    }
}

/**
 *  Move Ira along a wline
 */
void Markov::MoveWormOnW()
{
    if (!Worm->Exist)
        return;
    
    vertex Ira = Worm->Ira;
    vertex Masha = Worm->Masha;
    
    vertex v1 = Ira;
    wLine w = v1->NeighW();
    vertex v2 = w->NeighVer(INVERSE(v1->Dir));
    if (v2 == Ira || v2 == Masha)
        return;
    Momentum k = w->K + SIGN(v1->Dir) * Worm->K;
    //TODO: Check Hash for k

    Complex wWeight = W->Weight(v1->Dir, v1->R, v2->R, v1->Tau, v2->Tau, v1->Spin(),
                                v2->Spin(), w->IsWorm, w->IsMeasure, w->IsDelta);

    Complex weightRatio = wWeight / w->Weight;
    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    real wormWeight = WormWeight->Weight(v2->R, Masha->R, v2->Tau, Masha->Tau);
    prob *= wormWeight / Worm->Weight;

    if (prob >= 1.0 || RNG->urn() < prob) {
        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;

        w->Weight = wWeight;
        w->K = k;

        Ira = v2;
        Worm->Weight = wormWeight;
    }
}

/**
 *  reconnect from I->A, M->B to I->B, M->A
 */
void Markov::Reconnect()
{
    if (!Worm->Exist)
        return;
    
    vertex Ira = Worm->Ira;
    vertex Masha = Worm->Masha;
    
    if (!(Ira->R == Masha->R))
        return;
    int dir = RandomPickDir();
    gLine GIA = Ira->NeighG(dir);
    gLine GMB = Masha->NeighG(dir);
    if (GIA->Spin() != GMB->Spin())
        return;
    
    Momentum k = Worm->K + SIGN(dir) * (GMB->K - GIA->K);
    //TODO: Hash check for k

    vertex vA = GIA->NeighVer(dir);
    Complex GIAWeight = G->Weight(dir, Masha->R, vA->R, Masha->Tau, vA->Tau,
                                  Masha->Spin(dir), vA->Spin(INVERSE(dir)), GIA->IsMeasure);

    vertex vB = GMB->NeighVer(dir);
    Complex GMBWeight = G->Weight(dir, Ira->R, vB->R, Ira->Tau, vB->Tau,
                                  Ira->Spin(dir), vB->Spin(INVERSE(dir)), GMB->IsMeasure);

    Complex weightRatio = (-1) * GIAWeight * GMBWeight / (GIA->Weight * GMB->Weight);
    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);
    if (prob >= 1.0 || RNG->urn() < prob) {
        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;
        Diag->SignFermiLoop *= -1;

        Masha->nG[dir] = GIA;
        GIA->nVer[INVERSE(dir)] = Masha;

        Ira->nG[dir] = GMB;
        GMB->nVer[INVERSE(dir)] = Ira;

        Worm->K = k;

        GIA->Weight = GIAWeight;
        GMB->Weight = GMBWeight;
    }
}

void Markov::AddInteraction()
{
    if(!Worm->Exist)
        return;
    if(Diag->Order >= Order)
        return;
    vertex Ira = Worm->Ira, Masha = Worm->Masha;
    
    Momentum kW = RandomPickK();
    //TODO Hash Check for k
    
    int dir = RandomPickDir();
    int dirW = RandomPickDir();
    gLine GIC = Ira->NeighG(dir), GMD = Masha->NeighG(dir);
    
    Momentum kIA = GIC->K - SIGN(dir)*SIGN(dirW)*kW;
    Momentum kMB = GMD->K + SIGN(dir)*SIGN(dirW)*kW;
    Momentum kWorm = Worm->K - SIGN(dirW) * kW;
    //TODO: Hash Check for KIA, kMB, kWorm

    bool isdelta = RandomPickBool();
    real tauA=RandomPickTau(), tauB;
    
    if(isdelta)
        tauB = tauA;
    else
        tauB = RandomPickTau();
    
    spin spinA[2] = {GIC->Spin(), GIC->Spin()};
    spin spinB[2] = {GMD->Spin(), GMD->Spin()};
    vertex vC = GIC->NeighVer(dir), vD = GMD->NeighVer(dir);
    Site RA = vC->R, RB = vD->R;
    
    Complex wWeight = W->Weight(dirW, RA, RB, tauA, tauB, spinA, spinB, false, false, isdelta);
    Complex GIAWeight = G->Weight(FLIP(dir), Ira->R, RA, Ira->Tau, tauA,
                                  Ira->Spin(dir), spinA[FLIP(dir)], false);
    Complex GMBWeight = G->Weight(FLIP(dir), Masha->R, RB, Masha->Tau, tauB,
                                  Masha->Spin(dir), spinB[FLIP(dir)], false);
    Complex GACWeight = G->Weight(FLIP(dir), RA, vC->R, tauA, vC->Tau,
                                  spinA[dir], vC->Spin(FLIP(dir)), GIC->IsMeasure);
    Complex GBDWeight = G->Weight(FLIP(dir), RB, vD->R, tauB, vD->Tau,
                                  spinB[dir], vD->Spin(FLIP(dir)), GMD->IsMeasure);
    
    Complex weightRatio = (-1) * GIAWeight * GMBWeight *wWeight * GACWeight
                          *GBDWeight/(GIC->Weight * GMD->Weight);
    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);
    if (prob >= 1.0 || RNG->urn() < prob) {
        Diag->Order += 1;
        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;
        
        vertex vA = Diag->Ver.Add();
        vertex vB = Diag->Ver.Add();
        gLine GIA = Diag->G.Add();
        gLine GMB = Diag->G.Add();
        wLine WAB = Diag->W.Add();
        
        gLine ng[2];
        vertex nver[2];
        
        ng[dir] = GIC;
        ng[FLIP(dir)] = GIA;
        vA->SetVertex(RA, tauA, spinA, dirW, ng, WAB);
        
        ng[dir] = GMD;
        ng[FLIP(dir)] = GMB;
        vB->SetVertex(RB, tauB, spinB, FLIP(dirW), ng, WAB);
        
        nver[FLIP(dir)] = Ira;
        nver[dir] = vA;
        GIA->SetGLine(kIA, GIAWeight, false, nver);
        
        nver[FLIP(dir)] = Masha;
        nver[dir] = vB;
        GMB->SetGLine(kMB, GMBWeight, false, nver);
        
        nver[dirW] = vA;
        nver[FLIP(dirW)] = vB;
        WAB->SetWLine(kW, wWeight, false, false, isdelta, nver);
        
        Ira->nG[dir] = GIA;
        Masha->nG[dir] = GMB;
        GIC->nVer[FLIP(dir)] = vA;
        GMD->nVer[FLIP(dir)] = vB;
        
        Worm->K = kWorm;
        
        GIC->Weight = GACWeight;
        GMD->Weight = GBDWeight;
    }
}

void Markov::DeleteInteraction()
{
    if(!Worm->Exist)
        return;
    if(Diag->Order <= 1)
        return;
}



Momentum Markov::RandomPickK()
{
    return (Momentum)(RNG->irn(-MAX_K, MAX_K - 1));
}

int Markov::RandomPickDeltaSpin()
{
    return RNG->irn(0, 1) * 2 - 1;
}

int Markov::RandomPickDir()
{
    return RNG->irn(0, 1);
}

real Markov::RandomPickTau()
{
    return RNG->urn()*Beta;
}

bool Markov::RandomPickBool()
{
    return (RNG->irn(0,1)==0? true:false);
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
