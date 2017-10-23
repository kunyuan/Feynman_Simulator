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
using namespace weight;

#define SIGN(x) ((x) == IN ? 1 : -1)
#define NAME(x) #x

/**
 * Change the diagram from a Sigma to a dSigma/dG G^2 GammaG
 */
void Markov::JumpToGammaG() {
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

//    Site r_u = RandomPickSite();
    Site r_u = vA->R;

    spin spin_u = UP;

    UExt->R = r_u;
    UExt->Tau = tau_u;
    UExt->Spin = spin_u;

    Complex gammaGWeight = G->Weight(vA->R, vB->R, vA->Tau, vB->Tau, spinA, spinB, false, true, UExt);

    Complex weightRatio = gammaGWeight / (gAB->Weight);

    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    prob *= ProbofCall[JUMP_FROM_GAMMAG_TO_G] * 2.0* Diag->Order /
            (ProbofCall[JUMP_TO_GAMMAG] * ProbTau(tau_u));

//    prob *= ProbofCall[JUMP_FROM_GAMMAG_TO_G] / ProbofCall[JUMP_TO_GAMMAG];

    Proposed[JUMP_TO_GAMMAG][Diag->Order] += 1.0;
    if (prob >= 1.0 || RNG->urn() < prob) {
        Accepted[JUMP_TO_GAMMAG][Diag->Order] += 1.0;

        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;

        Diag->HasGammaGW = 1;

        Diag->Vin = &gAB->nVer[IN];
        Diag->Vout = &gAB->nVer[OUT];

        gAB->IsGammaG = true;
        gAB->Weight = gammaGWeight;
    }
}

/**
 * Change the diagram from a dSigma/dG G^2 GammaG to a Sigma
 */
void Markov::JumpFromGammaGToG()  {
    if (Diag->Order == 0 || Worm->Exist || !Diag->MeasureGLine || Diag->HasGammaGW!=1)
        return;

    vertex v_A = *Diag->Vin;
    vertex v_B = *Diag->Vout;
    ExtPoint * v_u = UExt;

    gLine gAB = v_A->NeighG(OUT);

    if (v_A->R != v_B->R)
        return;

    if (v_u->Spin != UP)
        return;

    if (v_A->R != v_u->R)
        return;

    Complex GABWeight = G->Weight(v_A->R, v_B->R, v_A->Tau, v_B->Tau,
                                  gAB->Spin(), gAB->Spin(),
                                  false, false, v_u);

    Complex gammaGWeight = gAB->Weight;

    Complex weightRatio = GABWeight/gammaGWeight;

    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    prob *= (ProbofCall[JUMP_TO_GAMMAG] * ProbTau(v_u->Tau))
                / (ProbofCall[JUMP_FROM_GAMMAG_TO_G] * 2.0 * Diag->Order);

    prob *= (ProbofCall[JUMP_TO_GAMMAG] ) / (ProbofCall[JUMP_FROM_GAMMAG_TO_G] );

    Proposed[JUMP_FROM_GAMMAG_TO_G][Diag->Order] += 1.0;
    if (prob >= 1.0 || RNG->urn() < prob) {
        Accepted[JUMP_FROM_GAMMAG_TO_G][Diag->Order] += 1.0;

        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;

        Diag->HasGammaGW = 0;

        gAB->Weight = GABWeight;
        gAB->IsGammaG = false;
    }
}

void Markov::JumpToGammaW() {
    return;
}

void Markov::JumpFromGammaWToW() {
}

void Markov::AddTwoG() {
    if (Diag->Order == 0 || Worm->Exist || !Diag->MeasureGLine || Diag->HasGammaGW != 1)
        return;

    if (Diag->MeasureGammaGW != 0)
        return;

    gLine measureG = Diag->GMeasure;
    vertex vA = measureG->NeighVer(IN);
    vertex vB = measureG->NeighVer(OUT);

    if (vA->R != vB->R)
        return;

    spin spin_meas = measureG->Spin();
    Site r_meas = vA->R;

    real tau_C = RandomPickTau();
    real tau_D = RandomPickTau();

    Complex GACWeight = G->Weight(r_meas, r_meas, vA->Tau, tau_C, spin_meas, spin_meas,
                                  false, false, UExt);

    Complex GDBWeight = G->Weight(r_meas, r_meas, tau_D, vB->Tau, spin_meas, spin_meas,
                                  false, false, UExt);
    Complex measureGWeight = G->Weight(r_meas, r_meas, tau_C, tau_D, spin_meas, spin_meas,
                                       true, false, UExt);

    Complex weightRatio = GACWeight * GDBWeight * measureGWeight/measureG->Weight;

    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    prob *= (ProbofCall[DELETE_TWO_G])
            / (ProbofCall[ADD_TWO_G] * ProbTau(tau_C) * ProbTau(tau_D));

//    real tau_C = vA->Tau;
//    real tau_D = vB->Tau;
//
//    Complex GACWeight = Complex(1.0, 0.0);
//    Complex GDBWeight = Complex(1.0, 0.0);

    ///Real G
//    Complex GACWeight = G->Weight(r_meas, r_meas, vA->Tau, tau_C, spin_meas, spin_meas,
//                                  false, false, UExt);
//
//    Complex GDBWeight = G->Weight(r_meas, r_meas, tau_D, vB->Tau, spin_meas, spin_meas,
//                                  false, false, UExt);

//    Complex measureGWeight = G->Weight(r_meas, r_meas, tau_C, tau_D, spin_meas, spin_meas,
//                                       true, false, UExt);
//
//    Complex weightRatio = GACWeight * GDBWeight * measureGWeight/measureG->Weight;
//
//    real prob = mod(weightRatio);
//    Complex sgn = phase(weightRatio);
//
//    prob *= (ProbofCall[DELETE_TWO_G]) / (ProbofCall[ADD_TWO_G]);

    Proposed[ADD_TWO_G][Diag->Order] += 1.0;
    if (prob >= 1.0 || RNG->urn() < prob) {
        Accepted[ADD_TWO_G][Diag->Order] += 1.0;

        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;

        Diag->MeasureGammaGW = 1;

        vertex vC = Diag->Ver.Add();
        vertex vD = Diag->Ver.Add();
        gLine GAC = Diag->G.Add();
        gLine GDB = Diag->G.Add();

        vC->nG[IN] = GAC;
        vC->nG[OUT] = measureG;
        spin vC_spin[2] = {spin_meas, spin_meas};
        vC->SetVertex(r_meas, tau_C, vC_spin, IN);

        vD->nG[IN] = measureG;
        vD->nG[OUT] = GDB;
        spin vD_spin[2] = {spin_meas, spin_meas};
        vD->SetVertex(r_meas, tau_D, vD_spin, OUT);

        GAC->nVer[IN] = vA;
        GAC->nVer[OUT] = vC;
        GAC->SetGLine(measureG->K, GACWeight,
                      false, false); //IsMeasure, IsGammaG

        GDB->nVer[IN] = vD;
        GDB->nVer[OUT] = vB;
        GDB->SetGLine(measureG->K, GDBWeight,
                      false, false); //IsMeasure, IsGammaG

        vA->nG[OUT] = GAC;
        vB->nG[IN] = GDB;
        measureG->nVer[IN] = vC;
        measureG->nVer[OUT] = vD;
        measureG->Weight = measureGWeight;

    }

}

void Markov::DeleteTwoG() {
    if (Diag->Order == 0 || Worm->Exist || !Diag->MeasureGLine || Diag->HasGammaGW != 1)
        return;

    if (Diag->MeasureGammaGW != 1)
        return;

    gLine measureG = Diag->GMeasure;

    vertex vC = measureG->NeighVer(IN);
    vertex vD = measureG->NeighVer(OUT);

    gLine GAC = vC->NeighG(IN);
    gLine GDB = vD->NeighG(OUT);

    vertex vA = GAC->NeighVer(IN);
    vertex vB = GDB->NeighVer(OUT);

    if (vC->R != vD->R || vC->R != vA->R || vD->R != vB->R)
        return;

    spin spin_meas = measureG->Spin();
    Site r_meas = vC->R;

    Complex measureGWeight = G->Weight(r_meas, r_meas, vA->Tau, vB->Tau, spin_meas, spin_meas, true, false, UExt);

    Complex weightRatio = measureGWeight/(measureG->Weight * GAC->Weight * GDB->Weight);

    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    prob *= (ProbofCall[ADD_TWO_G] * ProbTau(vC->Tau) * ProbTau(vD->Tau))/(ProbofCall[DELETE_TWO_G]);

//    if (vA->Tau != vC->Tau || vB->Tau != vD->Tau)
//        return;
//
//    Complex measureGWeight = G->Weight(r_meas, r_meas, vA->Tau, vB->Tau, spin_meas, spin_meas, true, false, UExt);
//
//    Complex weightRatio = measureGWeight/(measureG->Weight * GAC->Weight * GDB->Weight);
//
//    real prob = mod(weightRatio);
//    Complex sgn = phase(weightRatio);
//
//    prob *= (ProbofCall[ADD_TWO_G] )/(ProbofCall[DELETE_TWO_G]);

    Proposed[DELETE_TWO_G][Diag->Order] += 1.0;
    if (prob >= 1.0 || RNG->urn() < prob) {
        Accepted[DELETE_TWO_G][Diag->Order] += 1.0;

        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;

        Diag->MeasureGammaGW = 0;

        Diag->Ver.Remove(vC);
        Diag->Ver.Remove(vD);
        Diag->G.Remove(GAC);
        Diag->G.Remove(GDB);

        vA->nG[OUT] = measureG;
        vB->nG[IN] = measureG;
        measureG->nVer[IN] = vA;
        measureG->nVer[OUT] = vB;
        measureG->Weight = measureGWeight;
    }
}

void Markov::AddTwoW() {

}

void Markov::DeleteTwoW() {

}
