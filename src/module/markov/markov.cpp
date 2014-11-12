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
#include "module/weight/weight.h"
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
        else if (x < SumofProbofCall[7])
            ChangeTauOnVertex();
        else if (x < SumofProbofCall[8])
            ChangeROnVertex();
        else if (x < SumofProbofCall[9])
            ChangeRLoop();
        else if (x < SumofProbofCall[10])
            ChangeMeasureFromGToW();
        else if (x < SumofProbofCall[11])
            ChangeMeasureFromWToG();
        else if (x < SumofProbofCall[12])
            ChangeDeltaToContinuous();
        else if (x < SumofProbofCall[13])
            ChangeContinuousToDelta();
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
                                vin->Spin(), vout->Spin(),
                                true, //IsWorm
                                w->IsMeasure, w->IsDelta);
    
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
    vertex &Ira = Worm->Ira;
    vertex &Masha = Worm->Masha;
    
    wLine w = Ira->NeighW();
    if (!(w == Masha->NeighW()))
        return;
    Momentum k = SIGN(Ira->Dir)* w->K + Worm->K;
    //TODO: Hash check for k

    Complex wWeight = W->Weight(Ira->Dir, Ira->R, Masha->R, Ira->Tau, Masha->Tau,
                                Ira->Spin(), Masha->Spin(),
                                false, //IsWorm
                                w->IsMeasure, w->IsDelta);

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
    
    vertex &Ira = Worm->Ira;
    vertex &Masha = Worm->Masha;
    
    int dir = RandomPickDir();
    gLine g = Ira->NeighG(dir);
    vertex v2 = g->NeighVer(dir);

    if (v2 == Ira || v2 == Masha || CanNotMoveWorm(Worm->dSpin, g->Spin(), dir))
        return;
    Momentum k = g->K - SIGN(dir) * Worm->K;
    //TODO: Hash check for G(knew)
    //    if(!IsKValidG(k)) return;

    wLine w1 = Ira->NeighW();
    vertex vW1 = w1->NeighVer(INVERSE(Ira->Dir));
    bool isWormW1;
    if (vW1 == v2)
        isWormW1 = true;
    else
        isWormW1 = Diag->IsWorm(vW1);

    spin spinV1[2] = {Ira->Spin(0), Ira->Spin(1)};
    spinV1[dir] = FLIP(spinV1[dir]);

    Complex w1Weight = W->Weight(Ira->Dir, Ira->R, vW1->R, Ira->Tau, vW1->Tau,
                                 spinV1, vW1->Spin(), isWormW1, w1->IsMeasure, w1->IsDelta);

    wLine w2 = v2->NeighW();
    vertex vW2 = w2->NeighVer(INVERSE(v2->Dir));

    spin spinV2[2] = {v2->Spin(0), v2->Spin(1)};
    spinV2[INVERSE(dir)] = FLIP(spinV2[INVERSE(dir)]);

    Complex w2Weight = W->Weight(v2->Dir, v2->R, vW2->R, v2->Tau, vW2->Tau,
                                 spinV2, vW2->Spin(),
                                 true,  //IsWorm
                                 w2->IsMeasure, w2->IsDelta);

    Complex gWeight = G->Weight(INVERSE(dir), Ira->R, v2->R, Ira->Tau, v2->Tau,
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

        Ira->SetSpin(spinV1);
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
    
    vertex &Ira = Worm->Ira;
    vertex &Masha = Worm->Masha;
    
    wLine w = Ira->NeighW();
    vertex v2 = w->NeighVer(INVERSE(Ira->Dir));
    if (v2 == Ira || v2 == Masha)
        return;
    Momentum k = w->K + SIGN(Ira->Dir) * Worm->K;
    //TODO: Check Hash for k

    Complex wWeight = W->Weight(Ira->Dir, Ira->R, v2->R, Ira->Tau, v2->Tau, Ira->Spin(),
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

/**
 *  add interaction
 *  I---C  ===>  I---A---C
 *  M---D  ===>  M---B---D
 */
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
    
    Complex wWeight = W->Weight(dirW, RA, RB, tauA, tauB, spinA, spinB,
                                false,      //IsWorm
                                false,      //IsMeasure
                                isdelta);
    
    Complex GIAWeight = G->Weight(INVERSE(dir), Ira->R, RA, Ira->Tau, tauA,
                                  Ira->Spin(dir), spinA[INVERSE(dir)],
                                  false);  //IsMeasure
    
    Complex GMBWeight = G->Weight(INVERSE(dir), Masha->R, RB, Masha->Tau, tauB,
                                  Masha->Spin(dir), spinB[INVERSE(dir)],
                                  false);    //IsMeasure
    
    Complex GACWeight = G->Weight(INVERSE(dir), RA, vC->R, tauA, vC->Tau,
                                  spinA[dir], vC->Spin(INVERSE(dir)), GIC->IsMeasure);
    
    Complex GBDWeight = G->Weight(INVERSE(dir), RB, vD->R, tauB, vD->Tau,
                                  spinB[dir], vD->Spin(INVERSE(dir)), GMD->IsMeasure);
    
    Complex weightRatio = (-1) * GIAWeight * GMBWeight *wWeight * GACWeight
                          *GBDWeight/(GIC->Weight * GMD->Weight);
    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);
    
    if(isdelta)
        prob *= OrderWeight[Diag->Order+1]*ProbofCall[6]
            /(ProbofCall[5]*OrderWeight[Diag->Order]*ProbTau(tauA));
    else
        prob *= OrderWeight[Diag->Order+1]*ProbofCall[6]
            /(ProbofCall[5]*OrderWeight[Diag->Order]*ProbTau(tauA)*ProbTau(tauB));
    
    if (prob >= 1.0 || RNG->urn() < prob) {
        Diag->Order += 1;
        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;
        
        vertex vA = Diag->Ver.Add();
        vertex vB = Diag->Ver.Add();
        gLine GIA = Diag->G.Add();
        gLine GMB = Diag->G.Add();
        wLine WAB = Diag->W.Add();
        
        vA->nG[dir] = GIC;
        vA->nG[INVERSE(dir)] = GIA;
        vA->nW = WAB;
        vA->SetVertex(RA, tauA, spinA, dirW);
        
        vB->nG[dir] = GMD;
        vB->nG[INVERSE(dir)] = GMB;
        vB->nW = WAB;
        vB->SetVertex(RB, tauB, spinB, INVERSE(dirW));
        
        GIA->nVer[INVERSE(dir)] = Ira;
        GIA->nVer[dir] = vA;
        GIA->SetGLine(kIA, GIAWeight,
                      false);   //IsMeasure
        
        GMB->nVer[INVERSE(dir)] = Masha;
        GMB->nVer[dir] = vB;
        GMB->SetGLine(kMB, GMBWeight,
                      false);   //IsMeasure
        
        WAB->nVer[dirW] = vA;
        WAB->nVer[INVERSE(dirW)] = vB;
        WAB->SetWLine(kW, wWeight,
                      false,   //IsWorm
                      false,   //IsMeasure
                      isdelta);
        
        Ira->nG[dir] = GIA;
        Masha->nG[dir] = GMB;
        GIC->nVer[INVERSE(dir)] = vA;
        GMD->nVer[INVERSE(dir)] = vB;
        
        Worm->K = kWorm;
        
        GIC->Weight = GACWeight;
        GMD->Weight = GBDWeight;
    }
}

/**
 *  delete interaction
 *  I---A---C  ===>  I---C
 *  M---B---D  ===>  M---D
 */
void Markov::DeleteInteraction()
{
    if(!Worm->Exist)
        return;
    if(Diag->Order <= 1)
        return;
    vertex Ira = Worm->Ira, Masha = Worm->Masha;
    
    int dir = RandomPickDir();
    gLine GIA = Ira->NeighG(dir), GMB = Masha->NeighG(dir);
    if(GIA->IsMeasure) return;
    if(GMB->IsMeasure) return;
    
    vertex vA = GIA->NeighVer(dir), vB = GMB->NeighVer(dir);
    if(vA->Spin(IN)!=vA->Spin(OUT)) return;
    if(vB->Spin(IN)!=vB->Spin(OUT)) return;
    
    if(vA->NeighW()!=vB->NeighW()) return;
    
    wLine wAB = vA->NeighW();
    if(wAB->IsMeasure) return;
    if(wAB->IsWorm) return;
    
    gLine GAC = vA->NeighG(dir), GBD = vB->NeighG(dir);
    vertex vC = GAC->NeighVer(dir), vD = GBD->NeighVer(dir);
    if(vA->R != vC->R) return;
    if(vB->R != vD->R) return;
    
    Momentum kWorm = Worm->K + SIGN(vA->Dir) * wAB->K;
    //TODO: Hash Check for kWorm
    
    Complex GICWeight = G->Weight(INVERSE(dir), Ira->R, vC->R, Ira->Tau, vC->Tau,
                                  Ira->Spin(dir), vC->Spin(INVERSE(dir)), GAC->IsMeasure);
    Complex GMDWeight = G->Weight(INVERSE(dir), Masha->R, vD->R, Masha->Tau, vD->Tau,
                                  Masha->Spin(dir), vD->Spin(INVERSE(dir)), GBD->IsMeasure);
    
    Complex weightRatio = (-1) * GICWeight * GMDWeight/(GIA->Weight * GMB->Weight *GAC->Weight
                                                        * GBD->Weight *wAB->Weight);
    
    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);
    
    if(wAB->IsDelta)
        prob *= OrderWeight[Diag->Order]*ProbofCall[5]*ProbTau(vA->Tau)
                /(ProbofCall[6]*OrderWeight[Diag->Order+1]);
    else
        prob *= OrderWeight[Diag->Order]*ProbofCall[5]*ProbTau(vA->Tau)*ProbTau(vB->Tau)
                /(ProbofCall[6]*OrderWeight[Diag->Order+1]);
    
    if (prob >= 1.0 || RNG->urn() < prob) {
        Diag->Order -= 1;
        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;
        
        Diag->Ver.Remove(vA);
        Diag->Ver.Remove(vB);
        Diag->G.Remove(GIA);
        Diag->G.Remove(GMB);
        Diag->W.Remove(wAB);
        
        Ira->nG[dir] = GAC;
        Masha->nG[dir] = GBD;
        GAC->nVer[INVERSE(dir)] = Ira;
        GBD->nVer[INVERSE(dir)] = Masha;
        
        Worm->K = kWorm;
        
        GAC->Weight = GICWeight;
        GBD->Weight = GMDWeight;
    }
}

/**
 *  change Tau in a vertex
 */
void Markov::ChangeTauOnVertex()
{
    if(Worm->Exist)
        return;
    vertex ver = Diag->RandomPickVer();
    real tau = RandomPickTau();
    
    gLine gin = ver->NeighG(IN), gout = ver->NeighG(OUT);
    Complex ginWeight = G->Weight(gin->NeighVer(IN)->R, ver->R,
                                  gin->NeighVer(IN)->Tau, tau,
                                  gin->NeighVer(IN)->Spin(OUT), ver->Spin(IN),
                                  gin->IsMeasure);
    
    Complex goutWeight = G->Weight(ver->R, gout->NeighVer(OUT)->R,
                                  tau, gout->NeighVer(OUT)->Tau,
                                  ver->Spin(OUT), gout->NeighVer(OUT)->Spin(IN),
                                  gout->IsMeasure);
    
    wLine w = ver->NeighW();
    vertex vW = w->NeighVer(INVERSE(ver->Dir));
    Complex wWeight;
    if(w->IsDelta)
        wWeight = W->Weight(ver->Dir, ver->R, vW->R, tau, tau, ver->Spin(), vW->Spin(),
                            w->IsWorm, w->IsMeasure, w->IsDelta);
    else
        wWeight = W->Weight(ver->Dir, ver->R, vW->R, tau, vW->Tau, ver->Spin(), vW->Spin(),
                            w->IsWorm, w->IsMeasure, w->IsDelta);
    
    Complex weightRatio = ginWeight * goutWeight *wWeight
                        /(gin->Weight * gout->Weight * w->Weight);
    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);
    
    prob *= ProbTau(ver->Tau)/ProbTau(tau);
    
    if (prob >= 1.0 || RNG->urn() < prob) {
        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;
        
        ver->Tau = tau;
        if(w->IsDelta) vW->Tau = tau;
        
        gin->Weight = ginWeight;
        gout->Weight = goutWeight;
        w->Weight = wWeight;
    }
}

/**
 *  change R in a vertex
 */
void Markov::ChangeROnVertex()
{
    if(Worm->Exist)
        return;
    vertex ver = Diag->RandomPickVer();
    Site site = RandomPickSite();
    gLine gin = ver->NeighG(IN), gout = ver->NeighG(OUT);
    Complex ginWeight = G->Weight(gin->NeighVer(IN)->R, site,
                                  gin->NeighVer(IN)->Tau, ver->Tau,
                                  gin->NeighVer(IN)->Spin(OUT), ver->Spin(IN),
                                  gin->IsMeasure);
    
    Complex goutWeight = G->Weight(site, gout->NeighVer(OUT)->R,
                                  ver->Tau, gout->NeighVer(OUT)->Tau,
                                  ver->Spin(OUT), gout->NeighVer(OUT)->Spin(IN),
                                  gout->IsMeasure);
    wLine w = ver->NeighW();
    vertex vW = w->NeighVer(INVERSE(ver->Dir));
    Complex wWeight = W->Weight(ver->Dir, site, vW->R, ver->Tau, vW->Tau, ver->Spin(), vW->Spin(),
                            w->IsWorm, w->IsMeasure, w->IsDelta);
    
    Complex weightRatio = ginWeight * goutWeight *wWeight
                        /(gin->Weight * gout->Weight * w->Weight);
    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);
    
    prob *= ProbSite(ver->R)/ProbSite(site);
    
    if (prob >= 1.0 || RNG->urn() < prob) {
        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;
        
        ver->R = site;
        
        gin->Weight = ginWeight;
        gout->Weight = goutWeight;
        w->Weight = wWeight;
    }
}


/**
 *  change all R in a fermi loop
 */
void Markov::ChangeRLoop()
{
    if(Worm->Exist)
        return;
    //TODO: If G is not a local function, return;
    
    vertex v[2*MAX_ORDER]={nullptr};
    bool flagVer[2*MAX_ORDER]={false};
    int flagW[MAX_ORDER]={0};
    int n = 0;
    v[0] = Diag->RandomPickVer();
    Site oldR = v[0]->R;
    
    while(!flagVer[n])
    {
        flagVer[n]=true;
        flagW[v[n]->NeighW()->Name] ++;
        v[n+1] = v[n]->NeighG(OUT)->NeighVer(OUT);
        if(v[n+1]->R != oldR)
            return;
        n++;
    }
    
    Site newR = RandomPickSite();
    
    gLine g=nullptr;
    wLine w=nullptr;
    
    Complex GWeight[2*MAX_ORDER]={Complex(1.0, 0.0)};
    Complex WWeight[2*MAX_ORDER]={Complex(1.0, 0.0)};
    
    Complex oldWeight(1.0, 0.0);
    Complex newWeight(1.0, 0.0);
    for(int i=0; i<n; i++)
    {
        g = v[i]->NeighG(OUT);
        GWeight[i] = G->Weight(newR, newR, v[i]->Tau, g->NeighVer(OUT)->Tau,
                               v[i]->Spin(OUT), g->NeighVer(OUT)->Spin(IN), g->IsMeasure);
        newWeight *= GWeight[i];
        oldWeight *= g->Weight;
        
        w = v[i]->NeighW();
        if(flagW[w->Name]==1)
        {
            WWeight[i] *= W->Weight(v[i]->Dir, newR, w->NeighVer(INVERSE(v[i]->Dir))->R,
                                   v[i]->Tau, w->NeighVer(INVERSE(v[i]->Dir))->Tau,
                                   v[i]->Spin(), w->NeighVer(INVERSE(v[i]->Dir))->Spin(),
                                   w->IsWorm, w->IsMeasure, w->IsDelta);
            newWeight *= WWeight[i];
            oldWeight *= w->Weight;
        }else if(flagW[w->Name]==2){
            flagW[w->Name] = 0;
            WWeight[i] *= W->Weight(v[i]->Dir, newR, newR,
                                   v[i]->Tau, w->NeighVer(INVERSE(v[i]->Dir))->Tau,
                                   v[i]->Spin(), w->NeighVer(INVERSE(v[i]->Dir))->Spin(),
                                   w->IsWorm, w->IsMeasure, w->IsDelta);
            newWeight *= WWeight[i];
            oldWeight *= w->Weight;
        }
    }
    
    Complex weightRatio = newWeight/oldWeight;
    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);
    
    prob *= ProbSite(oldR)/ProbSite(newR);
    
    if (prob >= 1.0 || RNG->urn() < prob) {
        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;
        
        for(int i=0; i<n; i++)
        {
            v[i]->R = newR;
            v[i]->NeighG(OUT)->Weight = GWeight[i];
            v[i]->NeighW()->Weight = WWeight[i];
        }
    }
}


/**
 *  change measuring line from G to W
 */
void Markov::ChangeMeasureFromGToW()
{
    if(!Diag->MeasureGLine) return;
    
    wLine w = Diag->RandomPickW();
    if(w->IsDelta) return;
    
    gLine g = Diag->GMeasure;
    Complex gWeight = G->Weight(g->NeighVer(IN)->R, g->NeighVer(OUT)->R,
                                g->NeighVer(IN)->Tau, g->NeighVer(OUT)->Tau,
                                g->Spin(), g->Spin(),
                                false);  //IsMeasure
    
    Complex wWeight = W->Weight(w->NeighVer(IN)->R, w->NeighVer(OUT)->R,
                                w->NeighVer(IN)->Tau, w->NeighVer(OUT)->Tau,
                                w->NeighVer(IN)->Spin(), w->NeighVer(OUT)->Spin(),
                                w->IsWorm,
                                true,    //IsMeasure
                                w->IsDelta);
    
    Complex weightRatio = gWeight * wWeight
                        /(g->Weight * w->Weight);
    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);
    
    prob *= 0.50*ProbofCall[11]/(ProbofCall[10]);
    
    if (prob >= 1.0 || RNG->urn() < prob) {
        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;
        
        Diag->MeasureGLine = false;
        Diag->GMeasure = nullptr;
        Diag->WMeasure = w;
        
        g->IsMeasure = false;
        w->IsMeasure = true;
        
        g->Weight = gWeight;
        w->Weight = wWeight;
    }
}

/**
 *  change measuring line from W to G
 */
void Markov::ChangeMeasureFromWToG()
{
    if(Diag->MeasureGLine) return;
    
    gLine g = Diag->RandomPickG();
    
    wLine w = Diag->WMeasure;
    Complex gWeight = G->Weight(g->NeighVer(IN)->R, g->NeighVer(OUT)->R,
                                g->NeighVer(IN)->Tau, g->NeighVer(OUT)->Tau,
                                g->Spin(), g->Spin(),
                                true); //IsMeasure
    
    Complex wWeight = W->Weight(w->NeighVer(IN)->R, w->NeighVer(OUT)->R,
                                w->NeighVer(IN)->Tau, w->NeighVer(OUT)->Tau,
                                w->NeighVer(IN)->Spin(), w->NeighVer(OUT)->Spin(),
                                w->IsWorm,
                                false,   //IsMeasure
                                w->IsDelta);
    
    Complex weightRatio = gWeight * wWeight
                        /(g->Weight * w->Weight);
    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);
    
    prob *= ProbofCall[10]*2.0/(ProbofCall[11]);
    
    if (prob >= 1.0 || RNG->urn() < prob) {
        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;
        
        Diag->MeasureGLine = true;
        Diag->GMeasure = g;
        Diag->WMeasure = nullptr;
        
        g->IsMeasure = true;
        w->IsMeasure = false;
        
        g->Weight = gWeight;
        w->Weight = wWeight;
    }
}


/**
 *  change a Wline from Delta to continuous
 */
void Markov::ChangeDeltaToContinuous()
{
    wLine w = Diag->RandomPickW();
    vertex vin = w->NeighVer(IN), vout = w->NeighVer(OUT);
    gLine G1 = vout->NeighG(IN), G2 = vout->NeighG(OUT);
    if(!w->IsDelta) return;
    real tau = RandomPickTau();
    Complex wWeight = W->Weight(vin->R, vout->R, vin->Tau, tau, vin->Spin(),
                        vout->Spin(), w->IsWorm, w->IsMeasure,
                                false); //IsDelta
    
    Complex G1Weight = G->Weight(G1->NeighVer(IN)->R, vout->R,
                                 G1->NeighVer(IN)->Tau, tau,
                                 G1->NeighVer(IN)->Spin(OUT), vout->Spin(IN),
                                 G1->IsMeasure);
    
    Complex G2Weight = G->Weight(OUT, G2->NeighVer(OUT)->R, vout->R,
                                 G2->NeighVer(OUT)->Tau, tau,
                                 G2->NeighVer(OUT)->Spin(IN), vout->Spin(OUT),
                                 G2->IsMeasure);
    
    Complex weightRatio = G1Weight * G2Weight *wWeight
                        /(G1->Weight * G2->Weight * w->Weight);
    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    prob *= ProbofCall[13]/(ProbofCall[12]*ProbTau(tau));
    
    if (prob >= 1.0 || RNG->urn() < prob) {
        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;
        
        vout->Tau = tau;
        w->IsDelta = false;
        
        G1->Weight = G1Weight;
        G2->Weight = G2Weight;
        w->Weight = wWeight;
    }
}


/**
 *  Change a Wline from continuous to Delta
 */
void Markov::ChangeContinuousToDelta()
{
    wLine w = Diag->RandomPickW();
    if(w->IsDelta||w->IsMeasure) return;
    
    vertex vin = w->NeighVer(IN), vout = w->NeighVer(OUT);
    gLine G1 = vout->NeighG(IN), G2 = vout->NeighG(OUT);
    
    Complex wWeight = W->Weight(vin->R, vout->R, vin->Tau, vin->Tau, vin->Spin(),
                        vout->Spin(), w->IsWorm, w->IsMeasure,
                                true);   //IsDelta
    
    Complex G1Weight = G->Weight(G1->NeighVer(IN)->R, vout->R,
                                 G1->NeighVer(IN)->Tau, vin->Tau,
                                 G1->NeighVer(IN)->Spin(OUT), vout->Spin(IN),
                                 G1->IsMeasure);
    
    Complex G2Weight = G->Weight(OUT, G2->NeighVer(OUT)->R, vout->R,
                                 G2->NeighVer(OUT)->Tau, vin->Tau,
                                 G2->NeighVer(OUT)->Spin(IN), vout->Spin(OUT),
                                 G2->IsMeasure);
    
    Complex weightRatio = G1Weight * G2Weight *wWeight
                        /(G1->Weight * G2->Weight * w->Weight);
    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    prob *= ProbofCall[12]*ProbTau(vout->Tau)/ProbofCall[13];

    if (prob >= 1.0 || RNG->urn() < prob) {
        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;
        
        vout->Tau = vin->Tau;
        w->IsDelta = true;
        
        G1->Weight = G1Weight;
        G2->Weight = G2Weight;
        w->Weight = wWeight;
    }
}

Momentum Markov::RandomPickK()
{
    return (Momentum)(RNG->irn(-MAX_K + 1, MAX_K - 1));
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

real Markov::ProbTau(real tau)
{
    return 1.0/Beta;
}

bool Markov::RandomPickBool()
{
    return (RNG->irn(0,1)==0? true:false);
}

Site Markov::RandomPickSite()
{
    Vec<int> coord;
    for(int i=0; i<D; i++)
        coord[i] = RNG->irn(0, Lat->Size[i]-1);
    return (Site(RNG->irn(0,NSublattice-1), coord));
}
    
real Markov::ProbSite(const Site& site)
{
    
    return 1.0/(Lat->Vol*Lat->SublatVol);
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
