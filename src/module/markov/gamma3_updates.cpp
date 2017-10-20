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
#include "module/weight/gamma3.h"

using namespace std;
using namespace diag;
using namespace para;
using namespace mc;
using namespace weight;

#define SIGN(x) ((x) == IN ? 1 : -1)
#define NAME(x) #x

/**
 * Change the diagram from a Sigma to a dSigma/dG G^2 GammaG
 */
void Markov::JumpTodSdG() {
    if (Diag->MeasureGammaGW != 0)
        return;

    if (Diag->Order == 0 || Worm->Exist || !Diag->MeasureGLine || Diag->HasGammaGW != 0)
        return;

    gLine gAB = Diag->G.RandomPick(*RNG);

    if(gAB->IsMeasure)
        return;

    vertex vA = gAB->NeighVer(IN);
    vertex vB = gAB->NeighVer(OUT);

    spin spinA = gAB->Spin();
    spin spinB = gAB->Spin();

    real tau_u = RandomPickTau();
    Site r_u = RandomPickSite();
    spin spinu_out = RandomPickSpin();

    ExtPoint v_u = ExtPoint();
    v_u.R = r_u;
    v_u.Tau = tau_u;
    v_u.Spin = spinu_out;


    Complex gammaGWeight = G->Weight(vA->R, vB->R, vA->Tau, vB->Tau, spinA, spinB, false, true, &v_u);
//    Complex gammaGWeight = Complex(1.0, 0.0);

    Complex weightRatio = gammaGWeight / (gAB->Weight);

    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    prob *= ProbofCall[JUMP_FROM_DSDG_TO_SIGMA] * 2.0* Diag->Order / (ProbofCall[JUMP_TO_DSDG]
             * ProbTau(tau_u) * ProbSite(r_u) * 0.5);

    Proposed[JUMP_TO_DSDG][Diag->Order] += 1.0;
    if (prob >= 1.0 || RNG->urn() < prob) {
        Accepted[JUMP_TO_DSDG][Diag->Order] += 1.0;

        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;
        Diag->HasGammaGW = 1;

        Diag->Vin = &vA;
        Diag->Vout = &vB;
        Diag->V_Ext = v_u;

        gAB->IsGammaG = true;
        gAB->Weight = gammaGWeight;
    }
}

/**
 * Change the diagram from a dSigma/dG G^2 GammaG to a Sigma
 */
void Markov::JumpFromdSdGToSigma() {
    if (Diag->MeasureGammaGW != 0)
        return;
    if (Diag->Order == 0 || Worm->Exist || !Diag->MeasureGLine || Diag->HasGammaGW!=1)
        return;

    vertex v_A = *Diag->Vin;
    vertex v_B = *Diag->Vout;
    ExtPoint v_u = Diag->V_Ext;

    gLine gAB = v_A->NeighG(OUT);

    if (v_A->R != v_B->R)
        return;

    Complex GABWeight = G->Weight(v_A->R, v_B->R, v_A->Tau, v_B->Tau,
                                  gAB->Spin(), gAB->Spin(),
                                  false, false);

    Complex gammaGWeight = gAB->Weight;

//    Complex gammaGWeight = Complex(1.0, 0.0);

    Complex weightRatio = GABWeight/gammaGWeight;

    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    prob *= (ProbofCall[JUMP_TO_DSDG] * ProbTau(v_u.Tau) * ProbSite(v_u.R) * 0.5)
                / (ProbofCall[JUMP_FROM_DSDG_TO_SIGMA] * 2.0 * Diag->Order);

    Proposed[JUMP_FROM_DSDG_TO_SIGMA][Diag->Order] += 1.0;
    if (prob >= 1.0 || RNG->urn() < prob) {
        Accepted[JUMP_FROM_DSDG_TO_SIGMA][Diag->Order] += 1.0;

        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;
        Diag->HasGammaGW = 0;

        gAB->Weight = GABWeight;
        gAB->IsGammaG = false;
    }
}
