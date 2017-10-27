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
    if (Diag->Order == 0 || Worm->Exist || Diag->HasGammaGW != 0)
        return;

    gLine gAB = Diag->G.RandomPick(*RNG);

    if(gAB->IsMeasure)
        return;

    vertex vA = gAB->NeighVer(IN);
    vertex vB = gAB->NeighVer(OUT);

    spin spinA = gAB->Spin();
    spin spinB = gAB->Spin();

    real tau_u = RandomPickTau();

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
    if (Diag->Order == 0 || Worm->Exist || Diag->HasGammaGW!=1)
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
    if (Diag->Order == 0 || Worm->Exist || Diag->HasGammaGW != 0)
        return;

    wLine wAB = Diag->W.RandomPick(*RNG);

    if(wAB->IsMeasure)
        return;

    if(wAB->IsDelta)
        return;

    vertex vA = wAB->NeighVer(IN);
    vertex vB = wAB->NeighVer(OUT);

    real tau_u = RandomPickTau();
    Site r_u = RandomPickSite();
    spin spin_u = UP;
    UExt->R = r_u;
    UExt->Tau = tau_u;
    UExt->Spin = spin_u;

    Complex gammaWWeight = W->Weight(vA->R, vB->R, vA->Tau, vB->Tau, vA->Spin(), vB->Spin(),
                                     false, false, false, true, UExt);

    Complex weightRatio = gammaWWeight / (wAB->Weight);

    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    prob *= ProbofCall[JUMP_FROM_GAMMAW_TO_W] * Diag->Order /
            (ProbofCall[JUMP_TO_GAMMAW] * ProbTau(tau_u) * ProbSite(r_u));

    Proposed[JUMP_TO_GAMMAW][Diag->Order] += 1.0;
    if (prob >= 1.0 || RNG->urn() < prob) {
        Accepted[JUMP_TO_GAMMAW][Diag->Order] += 1.0;

        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;

        Diag->HasGammaGW = 2;

        Diag->Vin = &wAB->nVer[IN];
        Diag->Vout = &wAB->nVer[OUT];

        wAB->IsGammaW = true;
        wAB->Weight = gammaWWeight;
    }
}

void Markov::JumpFromGammaWToW() {
    if (Diag->Order == 0 || Worm->Exist || Diag->HasGammaGW!=2)
        return;

    vertex v_A = *Diag->Vin;
    vertex v_B = *Diag->Vout;
    ExtPoint * v_u = UExt;

    wLine wAB = v_A->NeighW();

    if (v_u->Spin != UP)
        return;

    Complex WABWeight = W->Weight(v_A->R, v_B->R, v_A->Tau, v_B->Tau,
                                  v_A->Spin(), v_B->Spin(),
                                  false, false, false, false, v_u);

    Complex gammaWWeight = wAB->Weight;

    Complex weightRatio = WABWeight/gammaWWeight;

    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    prob *= (ProbofCall[JUMP_TO_GAMMAW] * ProbTau(v_u->Tau) * ProbSite(v_u->R))
            / (ProbofCall[JUMP_FROM_GAMMAW_TO_W] * Diag->Order);

    Proposed[JUMP_FROM_GAMMAW_TO_W][Diag->Order] += 1.0;
    if (prob >= 1.0 || RNG->urn() < prob) {
        Accepted[JUMP_FROM_GAMMAW_TO_W][Diag->Order] += 1.0;

        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;

        Diag->HasGammaGW = 0;

        wAB->Weight = WABWeight;
        wAB->IsGammaW = false;
    }
}

void Markov::AddTwoG() {
    if (Diag->Order == 0 || Worm->Exist || !Diag->MeasureGLine || Diag->HasGammaGW == 0)
        return;

    if (Diag->MeasureGammaGW != 0)
        return;

    gLine measureG = Diag->GMeasure;
    vertex vA = measureG->NeighVer(IN);
    vertex vB = measureG->NeighVer(OUT);

    ASSERT_ALLWAYS(vA->R == vB->R, "Measuring G Line R doesn't match!");

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
    if (Diag->Order == 0 || Worm->Exist || !Diag->MeasureGLine || Diag->HasGammaGW == 0)
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
    ASSERT_ALLWAYS(vC->Spin()[IN]== spin_meas, "Spin doesn't match!");
    ASSERT_ALLWAYS(vD->Spin()[OUT]== spin_meas, "Spin doesn't match!");

    Site r_meas = vC->R;

    Complex measureGWeight = G->Weight(r_meas, r_meas, vA->Tau, vB->Tau, spin_meas, spin_meas, true, false, UExt);

    Complex weightRatio = measureGWeight/(measureG->Weight * GAC->Weight * GDB->Weight);

    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    prob *= (ProbofCall[ADD_TWO_G] * ProbTau(vC->Tau) * ProbTau(vD->Tau))/(ProbofCall[DELETE_TWO_G]);

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

    if (Diag->Order == 0 || Worm->Exist || Diag->MeasureGLine || Diag->HasGammaGW == 0)
        return;

    if (Diag->MeasureGammaGW != 0)
        return;

    wLine measureW = Diag->WMeasure;
    vertex vA = measureW->NeighVer(IN);
    vertex vB = measureW->NeighVer(OUT);

    bool flip_C = RandomPickBool();
    bool flip_D = RandomPickBool();

    spin spin_C[2], spin_D[2];

    if(flip_C){
        spin_C[IN] = FLIP(vA->Spin(IN));
        spin_C[OUT] = FLIP(vA->Spin(OUT));
    }else{
        spin_C[IN] = vA->Spin(IN);
        spin_C[OUT] = vA->Spin(OUT);
    }

    if(spin_C[IN]+vA->Spin(IN) != spin_C[OUT]+vA->Spin(OUT))
        return;

    if(flip_D){
        spin_D[IN] = FLIP(vB->Spin(IN));
        spin_D[OUT] = FLIP(vB->Spin(OUT));
    }else{
        spin_D[IN] = vB->Spin(IN);
        spin_D[OUT] = vB->Spin(OUT);
    }

    if(spin_D[IN]+vB->Spin(IN) != spin_D[OUT]+vB->Spin(OUT))
        return;

    Site r_C = RandomPickSite();
    Site r_D = RandomPickSite();

//    bool IsDelta_C = RandomPickBool();
//    bool IsDelta_D = RandomPickBool();

    //TODO Only Delta W
    bool IsDelta_C = true;
    bool IsDelta_D = true;

    //TODO Only Continuous W
//    bool IsDelta_C = false;
//    bool IsDelta_D = false;

    real tau_C, tau_D;
    Complex WACWeight, WDBWeight;

    if (IsDelta_C){
        tau_C = vA->Tau;
    }else{
        tau_C = RandomPickTau();
    }

    if (IsDelta_D){
        tau_D = vB->Tau;
    }else{
        tau_D = RandomPickTau();
    }

    WACWeight = W->Weight(vA->R, r_C, vA->Tau, tau_C, vA->Spin(), spin_C,
                          false, false, IsDelta_C, false, UExt);

    WDBWeight = W->Weight(r_D, vB->R, tau_D, vB->Tau, spin_D, vB->Spin(),
                          false, false, IsDelta_D, false, UExt);

    Complex measureWWeight = W->Weight(r_C, r_D, tau_C, tau_D, spin_C, spin_D,
                                       false, true, false, false, UExt);

    Complex weightRatio = WACWeight * WDBWeight * measureWWeight/measureW->Weight;

    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

//    prob *= (ProbofCall[DELETE_TWO_W])
//            / (ProbofCall[ADD_TWO_W] *0.5 * 0.5 *0.5 *0.5 * ProbSite(r_C) * ProbSite(r_D));

    //TODO: only delta/continuous W
    prob *= (ProbofCall[DELETE_TWO_W])
            / (ProbofCall[ADD_TWO_W] *0.5 *0.5 * ProbSite(r_C) * ProbSite(r_D));

    if (!IsDelta_C)
        prob /= ProbTau(tau_C);
    if (!IsDelta_D)
        prob /= ProbTau(tau_D);

    Proposed[ADD_TWO_W][Diag->Order] += 1.0;
    if (prob >= 1.0 || RNG->urn() < prob) {
        Accepted[ADD_TWO_W][Diag->Order] += 1.0;

        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;

        Diag->MeasureGammaGW = 2;

        vertex vC = Diag->Ver.Add();
        vertex vD = Diag->Ver.Add();
        wLine WAC = Diag->W.Add();
        wLine WDB = Diag->W.Add();

        vC->nW = WAC;
        vC->SetVertex(r_C, tau_C, spin_C, IN);

        vD->nW = WDB;
        vD->SetVertex(r_D, tau_D, spin_D, OUT);

        WAC->nVer[IN] = vA;
        WAC->nVer[OUT] = vC;
        WAC->SetWLine(measureW->K, WACWeight,
                      false, false, IsDelta_C, false); //IsMeasure, IsGammaG

        WDB->nVer[IN] = vD;
        WDB->nVer[OUT] = vB;
        WDB->SetWLine(measureW->K, WDBWeight,
                      false, false, IsDelta_D, false); //IsMeasure, IsGammaG

        vA->nW = WAC;
        vB->nW = WDB;
        measureW->nVer[IN] = vC;
        measureW->nVer[OUT] = vD;
        measureW->Weight = measureWWeight;

    }

}

void Markov::DeleteTwoW() {
    if (Diag->Order == 0 || Worm->Exist ||  Diag->MeasureGLine || Diag->HasGammaGW == 0)
        return;

    if (Diag->MeasureGammaGW != 2)
        return;

    wLine measureW = Diag->WMeasure;

    vertex vC = measureW->NeighVer(IN);
    vertex vD = measureW->NeighVer(OUT);

    wLine WAC = vC->NeighW();
    wLine WDB = vD->NeighW();

    vertex vA = WAC->NeighVer(IN);
    vertex vB = WDB->NeighVer(OUT);

    Complex measureWWeight = W->Weight(vA->R, vB->R, vA->Tau, vB->Tau, vA->Spin(), vB->Spin(), false, true, false, false, UExt);

    Complex weightRatio = measureWWeight/(measureW->Weight * WAC->Weight * WDB->Weight);

    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

//    prob *= (ProbofCall[ADD_TWO_W] * 0.5 *0.5 *0.5 *0.5 * ProbSite(vC->R) *ProbSite(vD->R))/(ProbofCall[DELETE_TWO_W]);

    //TODO only Delta/Continuous W
    prob *= (ProbofCall[ADD_TWO_W] *0.5 *0.5 * ProbSite(vC->R) *ProbSite(vD->R))/(ProbofCall[DELETE_TWO_W]);

    if (!WAC->IsDelta)
        prob *= ProbTau(vC->Tau);
    if (!WDB->IsDelta)
        prob *= ProbTau(vD->Tau);

    Proposed[DELETE_TWO_W][Diag->Order] += 1.0;
    if (prob >= 1.0 || RNG->urn() < prob) {
        Accepted[DELETE_TWO_W][Diag->Order] += 1.0;

        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;

        Diag->MeasureGammaGW = 0;

        Diag->Ver.Remove(vC);
        Diag->Ver.Remove(vD);
        Diag->W.Remove(WAC);
        Diag->W.Remove(WDB);

        vA->nW = measureW;
        vB->nW = measureW;
        measureW->nVer[IN] = vA;
        measureW->nVer[OUT] = vB;
        measureW->Weight = measureWWeight;
    }
}
