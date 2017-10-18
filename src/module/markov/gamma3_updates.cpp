//
// Created by yuanhuang on 10/18/17.
//
#include <stdio.h>
#include "markov.h"
#include "math.h"
#include "utility/utility.h"
#include "utility/momentum.h"
#include "module/diagram/diagram.h"
#include "module/parameter/parameter.h"
#include "lattice/lattice.h"
#include "module/weight/weight.h"
#include "module/weight/component.h"

using namespace std;
using namespace diag;
using namespace para;
using namespace mc;

#define SIGN(x) ((x) == IN ? 1 : -1)
#define NAME(x) #x

/**
 * Change the diagram from a Sigma to a dSigma/dG G^2 GammaG
 */
void Markov::JumpTodSdG() {
    if (Diag->Order == 0 || Worm->Exist || Diag->MeasureGLine)
        return;

    gLine gAB = Diag->G.RandomPick(*RNG);
    vertex vA = gAB->NeighVer(IN);
    vertex vB = gAB->NeighVer(OUT);

    real tauC = RandomPickTau(), tauD = RandomPickTau();

    spin spinAC = gAB->Spin();
    spin spinDB = gAB->Spin();
    Complex GACWeight = G->Weight(vA->R, vA->R, vA->Tau, tauC,
                                  spinAC, spinAC,
                                  false); //IsMeasure

    Complex GDBWeight = G->Weight(vB->R, vB->R, tauD, vB->Tau,
                                  spinDB, spinDB,
                                  false); //IsMeasure

    real tauE = RandomPickTau();
    Site rE = RandomPickSite();
    spin spinE = RandomPickSpin();
    Complex gammaG = GammaG->Weight(vA->R, vB->R, rE, tauC, tauD, tauE, spinAC, spinDB, spinE);

}


/**
 * Change the diagram from a dSigma/dG G^2 GammaG back to a Sigma
 */
void Markov::JumpFromdSdGToSigma() {
}


/**
 *  add delta interaction
 *  I---C  ===>  I---A---C
 *  M---D  ===>  M---B---D
 */
//void Markov::AddDeltaInteraction()
//{
//    if (!Worm->Exist)
//        return;
//    if (Diag->Order == 0 || Diag->Order >= Order)
//        return;
//    vertex Ira = Worm->Ira, Masha = Worm->Masha;
//
//    Momentum kW = RandomPickK();
//    if (Diag->WHashCheck(kW))
//        return;
//
//    int dir = RandomPickDir();
//    int dirW = RandomPickDir();
//    gLine GIC = Ira->NeighG(dir), GMD = Masha->NeighG(dir);
//
//    Momentum kIA = GIC->K - SIGN(dir) * SIGN(dirW) * kW;
//    Momentum kMB = GMD->K + SIGN(dir) * SIGN(dirW) * kW;
//    Momentum kWorm = Worm->K - SIGN(dirW) * kW;
//    if (Diag->GHashCheck(kIA))
//        return;
//    if (Diag->GHashCheck(kMB))
//        return;
//    if (kIA == kMB)
//        return;
//
//    real tauA = RandomPickTau();
//
//    spin spinA[2] = {GIC->Spin(), GIC->Spin()};
//    spin spinB[2] = {GMD->Spin(), GMD->Spin()};
//    vertex vC = GIC->NeighVer(dir), vD = GMD->NeighVer(dir);
//    Site RA = vC->R, RB = vD->R;
//
//    Complex wWeight = W->Weight(dirW, RA, RB, tauA, tauA, spinA, spinB,
//                                false, //IsWorm
//                                false, //IsMeasure
//                                true); //IsDelta
//
//    Complex GIAWeight = G->Weight(INVERSE(dir), Ira->R, RA, Ira->Tau, tauA,
//                                  Ira->Spin(dir), spinA[INVERSE(dir)],
//                                  false); //IsMeasure
//
//    Complex GMBWeight = G->Weight(INVERSE(dir), Masha->R, RB, Masha->Tau, tauA,
//                                  Masha->Spin(dir), spinB[INVERSE(dir)],
//                                  false); //IsMeasure
//
//    Complex GACWeight = G->Weight(INVERSE(dir), RA, vC->R, tauA, vC->Tau,
//                                  spinA[dir], vC->Spin(INVERSE(dir)), GIC->IsMeasure);
//
//    Complex GBDWeight = G->Weight(INVERSE(dir), RB, vD->R, tauA, vD->Tau,
//                                  spinB[dir], vD->Spin(INVERSE(dir)), GMD->IsMeasure);
//
//    Complex weightRatio = (-1) * GIAWeight * GMBWeight * wWeight * GACWeight * GBDWeight / (GIC->Weight * GMD->Weight);
//
//    real prob = mod(weightRatio);
//    Complex sgn = phase(weightRatio);
//
//    prob *= OrderReWeight[Diag->Order + 1] * ProbofCall[DEL_DELTA_INTERACTION] / (ProbofCall[ADD_DELTA_INTERACTION] * OrderReWeight[Diag->Order] * ProbTau(tauA));
//
//    Proposed[ADD_DELTA_INTERACTION][Diag->Order] += 1.0;
//    if (prob >= 1.0 || RNG->urn() < prob) {
//        Accepted[ADD_DELTA_INTERACTION][Diag->Order] += 1.0;
//        Diag->Order += 1;
//        Diag->Phase *= sgn;
//        Diag->Weight *= weightRatio;
//
//        vertex vA = Diag->Ver.Add();
//        vertex vB = Diag->Ver.Add();
//        gLine GIA = Diag->G.Add();
//        gLine GMB = Diag->G.Add();
//        wLine WAB = Diag->W.Add();
//
//        vA->nG[dir] = GIC;
//        vA->nG[INVERSE(dir)] = GIA;
//        vA->nW = WAB;
//        vA->SetVertex(RA, tauA, spinA, dirW);
//
//        vB->nG[dir] = GMD;
//        vB->nG[INVERSE(dir)] = GMB;
//        vB->nW = WAB;
//        vB->SetVertex(RB, tauA, spinB, INVERSE(dirW));
//
//        GIA->nVer[INVERSE(dir)] = Ira;
//        GIA->nVer[dir] = vA;
//        GIA->SetGLine(kIA, GIAWeight,
//                      false); //IsMeasure
//        Diag->AddGHash(kIA);
//
//        GMB->nVer[INVERSE(dir)] = Masha;
//        GMB->nVer[dir] = vB;
//        GMB->SetGLine(kMB, GMBWeight,
//                      false); //IsMeasure
//        Diag->AddGHash(kMB);
//
//        WAB->nVer[dirW] = vA;
//        WAB->nVer[INVERSE(dirW)] = vB;
//        WAB->SetWLine(kW, wWeight,
//                      false, //IsWorm
//                      false, //IsMeasure
//                      true); //IsDelta
//
//        Diag->AddWHash(kW);
//
//        Ira->nG[dir] = GIA;
//        Masha->nG[dir] = GMB;
//        GIC->nVer[INVERSE(dir)] = vA;
//        GMD->nVer[INVERSE(dir)] = vB;
//
//        Worm->K = kWorm;
//
//        GIC->Weight = GACWeight;
//        GMD->Weight = GBDWeight;
//    }
//}
