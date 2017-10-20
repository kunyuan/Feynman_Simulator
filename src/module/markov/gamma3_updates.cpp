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

#define SIGN(x) ((x) == IN ? 1 : -1)
#define NAME(x) #x

/**
 * Change the diagram from a Sigma to a dSigma/dG G^2 GammaG
 */
void Markov::JumpTodSdG() {
    if (Diag->Order == 0 || Worm->Exist || !Diag->MeasureGLine)
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

    Complex gammaGWeight = GammaG->Weight(vA->R, vB->R, r_u, vA->Tau, vB->Tau, tau_u, spinA, spinB, spinu_out);
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
        Diag->MeasuredSdG = true;

        vertex v_u = Diag->Ver.Add();
        Diag->Vin = &vA;
        Diag->Vout = &vB;
        Diag->Vw = &v_u;

        spin spin_u[2] = {spinu_out, spinu_out};
        v_u->SetVertex(r_u, tau_u, spin_u, OUT);
    }
}

/**
 * Change the diagram from a dSigma/dG G^2 GammaG to a Sigma
 */
void Markov::JumpFromdSdGToSigma() {
    if (Diag->Order == 0 || Worm->Exist || !Diag->MeasureGLine || !Diag->MeasuredSdG)
        return;

    vertex v_A = *Diag->Vin;
    vertex v_B = *Diag->Vout;
    vertex v_u = *Diag->Vw;

    gLine gAB = v_A->NeighG(OUT);

    if (v_A->R != v_B->R)
        return;

    Complex GABWeight = G->Weight(v_A->R, v_B->R, v_A->Tau, v_B->Tau,
                                  gAB->Spin(), gAB->Spin(),
                                  false); //IsMeasure

    Complex gammaGWeight = GammaG->Weight(v_A->R, v_B->R, v_u->R, v_A->Tau, v_B->Tau,
                                 v_u->Tau, gAB->Spin(), gAB->Spin(), v_u->NeighG(OUT)->Spin());
//    Complex gammaGWeight = Complex(1.0, 0.0);

    Complex weightRatio = GABWeight/gammaGWeight;

    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    prob *= (ProbofCall[JUMP_TO_DSDG] * ProbTau(v_u->Tau) * ProbSite(v_u->R) * 0.5)
                / (ProbofCall[JUMP_FROM_DSDG_TO_SIGMA] * 2.0 * Diag->Order);

    Proposed[JUMP_FROM_DSDG_TO_SIGMA][Diag->Order] += 1.0;
    if (prob >= 1.0 || RNG->urn() < prob) {
        Accepted[JUMP_FROM_DSDG_TO_SIGMA][Diag->Order] += 1.0;

        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;
        Diag->MeasuredSdG = false;

        gAB->Weight = GABWeight;

        Diag->Ver.Remove(v_u);
    }
}
