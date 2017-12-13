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
void Markov::JumpToGGGammaG() {
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

    Complex GGGammaGWeight = G->Weight(vA->R, vB->R, vA->Tau, vB->Tau, spinA, spinB, false, true, UExt);

    Complex weightRatio = GGGammaGWeight / (gAB->Weight);

    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    prob *= ProbofCall[JUMP_FROM_GGGAMMAG_TO_G] * 2.0* Diag->Order /
            (ProbofCall[JUMP_TO_GGGAMMAG] * ProbTau(tau_u));

    if(Diag->MeasureGLine)
        prob *= (*GammaGReweight);
    else
        prob *= (*GammaWReweight);


    Proposed[JUMP_TO_GGGAMMAG][Diag->Order] += 1.0;
    if (prob >= 1.0 || RNG->urn() < prob) {
        Accepted[JUMP_TO_GGGAMMAG][Diag->Order] += 1.0;

        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;

        Diag->HasGammaGW = 1;

        Diag->Vin = &gAB->nVer[IN];
        Diag->Vout = &gAB->nVer[OUT];

        gAB->IsGGGammaG = true;
        gAB->Weight = GGGammaGWeight;
    }
}

/**
 * Change the diagram from a dSigma/dG G^2 GammaG to a Sigma
 */
void Markov::JumpFromGGGammaGToG()  {
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

    Complex GGGammaGWeight = gAB->Weight;

    Complex weightRatio = GABWeight/GGGammaGWeight;

    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    prob *= (ProbofCall[JUMP_TO_GGGAMMAG] * ProbTau(v_u->Tau))
                / (ProbofCall[JUMP_FROM_GGGAMMAG_TO_G] * 2.0 * Diag->Order);

    if(Diag->MeasureGLine)
        prob /= (*GammaGReweight);
    else
        prob /= (*GammaWReweight);


    Proposed[JUMP_FROM_GGGAMMAG_TO_G][Diag->Order] += 1.0;
    if (prob >= 1.0 || RNG->urn() < prob) {
        Accepted[JUMP_FROM_GGGAMMAG_TO_G][Diag->Order] += 1.0;

        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;

        Diag->HasGammaGW = 0;

        gAB->Weight = GABWeight;
        gAB->IsGGGammaG = false;
    }
}

void Markov::JumpToWWGammaW() {
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

    Complex WWGammaWWeight = W->Weight(vA->R, vB->R, vA->Tau, vB->Tau, vA->Spin(), vB->Spin(),
                                     false, false, false, true, UExt);

    Complex weightRatio = WWGammaWWeight / (wAB->Weight);

    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    prob *= ProbofCall[JUMP_FROM_WWGAMMAW_TO_W] * Diag->Order /
            (ProbofCall[JUMP_TO_WWGAMMAW] * ProbTau(tau_u) * ProbSite(r_u));

    if(Diag->MeasureGLine)
        prob *= (*GammaGReweight);
    else
        prob *= (*GammaWReweight);

    Proposed[JUMP_TO_WWGAMMAW][Diag->Order] += 1.0;
    if (prob >= 1.0 || RNG->urn() < prob) {
        Accepted[JUMP_TO_WWGAMMAW][Diag->Order] += 1.0;

        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;

        Diag->HasGammaGW = 2;

        Diag->Vin = &wAB->nVer[IN];
        Diag->Vout = &wAB->nVer[OUT];

        wAB->IsWWGammaW = true;
        wAB->Weight = WWGammaWWeight;
    }
}

void Markov::JumpFromWWGammaWToW() {
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

    Complex WWGammaWWeight = wAB->Weight;

    Complex weightRatio = WABWeight/WWGammaWWeight;

    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    prob *= (ProbofCall[JUMP_TO_WWGAMMAW] * ProbTau(v_u->Tau) * ProbSite(v_u->R))
            / (ProbofCall[JUMP_FROM_WWGAMMAW_TO_W] * Diag->Order);

    if(Diag->MeasureGLine)
        prob /= (*GammaGReweight);
    else
        prob /= (*GammaWReweight);

    Proposed[JUMP_FROM_WWGAMMAW_TO_W][Diag->Order] += 1.0;
    if (prob >= 1.0 || RNG->urn() < prob) {
        Accepted[JUMP_FROM_WWGAMMAW_TO_W][Diag->Order] += 1.0;

        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;

        Diag->HasGammaGW = 0;

        wAB->Weight = WABWeight;
        wAB->IsWWGammaW = false;
    }
}

