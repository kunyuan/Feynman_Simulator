//
//  markov.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/17/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
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

bool CanNotMoveWorm(int dspin, spin sin, spin sout);
bool CanNotMoveWorm(int dspin, spin sin, int dir);

bool Markov::BuildNew(ParaMC &para, Diagram &diag, weight::Weight &weight)
{
    Reset(para, diag, weight);
    ASSERT_ALLWAYS(NUpdates >= (int)Operations::END,
                   "NUpdates " << NUpdates << " should larger than " << (int)Operations::END);

    InitialArray(ProbofCall, 0.0, NUpdates);
    InitialArray(SumofProbofCall, 0.0, NUpdates);

    for (int i = 0; i < NUpdates; i++) {
        ProbofCall[i] = 1.0 / real(NUpdates);
        for (int j = i; j < NUpdates; j++)
            SumofProbofCall[j] += ProbofCall[i];
    }

    InitialArray(&Accepted[0][0], 0.0, NUpdates * MAX_ORDER);
    InitialArray(&Proposed[0][0], 0.0, NUpdates * MAX_ORDER);

    OperationName[CREATE_WORM] = NAME(CREATE_WORM);
    OperationName[DELETE_WORM] = NAME(DELETE_WORM);
    OperationName[MOVE_WORM_G] = NAME(MOVE_WORM_G);
    OperationName[MOVE_WORM_W] = NAME(MOVE_WORM_W);
    OperationName[RECONNECT] = NAME(RECONNECT);
    OperationName[ADD_INTERACTION] = NAME(ADD_INTERACTION);
    OperationName[DEL_INTERACTION] = NAME(DEL_INTERACTION);
    OperationName[ADD_DELTA_INTERACTION] = NAME(ADD_DELTA_INTERACTION);
    OperationName[DEL_DELTA_INTERACTION] = NAME(DEL_DELTA_INTERACTION);
    OperationName[CHANGE_TAU_VERTEX] = NAME(CHANGE_TAU_VERTEX);
    OperationName[CHANGE_R_VERTEX] = NAME(CHANGE_R_VERTEX);
    OperationName[CHANGE_R_LOOP] = NAME(CHANGE_R_LOOP);
    OperationName[CHANGE_MEASURE_G2W] = NAME(CHANGE_MEASURE_G2W);
    OperationName[CHANGE_MEASURE_W2G] = NAME(CHANGE_MEASURE_W2G);
    OperationName[CHANGE_CONTINUS2DELTA] = NAME(CHANGE_CONTINUS2DELTA);
    OperationName[CHANGE_DELTA2CONTINUS] = NAME(CHANGE_DELTA2CONTINUS);
    OperationName[CHANGE_SPIN_VERTEX] = NAME(CHANGE_SPIN_VERTEX);
    OperationName[JUMP_TO_ORDER0] = NAME(JUMP_TO_ORDER0);
    OperationName[JUMP_BACK_TO_ORDER1] = NAME(JUMP_BACK_TO_ORDER1);
    return true;
}

void Markov::Reset(ParaMC &para, Diagram &diag, weight::Weight &weight)
{
    Beta = para.Beta;
    Order = para.Order;
    Counter = &para.Counter;
    Lat = &para.Lat;
    OrderReWeight = para.OrderReWeight.data();
    WormSpaceReweight = &para.WormSpaceReweight;
    PolarReweight = &para.PolarReweight;
    Diag = &diag;
    Worm = &diag.Worm;
    Sigma = weight.Sigma;
    Polar = weight.Polar;
    G = weight.G;
    W = weight.W;
    RNG = &para.RNG;
}

std::string Markov::_DetailBalanceStr(Operations op)
{
    string Output = OperationName[op] + ":\n";
    char temp[80];
    real TotalProposed = 0.0, TotalAccepted = 0.0;
    for (int i = 0; i <= Order; i++) {
        if (!Equal(Proposed[op][i], 0.0)) {
            TotalAccepted += Accepted[op][i];
            TotalProposed += Proposed[op][i];
            sprintf(temp, "\t%8s%2i:%15g%15g%15g\n", "Order", i, Proposed[op][i], Accepted[op][i], Accepted[op][i] / Proposed[op][i]);
            Output += temp;
        }
    }
    if (!Equal(TotalProposed, 0.0)) {
        sprintf(temp, "\t%10s:%15g%15g%15g\n", "Summation", TotalProposed, TotalAccepted, TotalAccepted / TotalProposed);
        Output += temp;
    }
    else
        Output += "\tNone\n";
    return Output;
}

std::string Markov::_CheckBalance(Operations op1, Operations op2)
{
    string Output = OperationName[op1] + "<------>" + OperationName[op2] + ":\n";
    char temp[80];
    real TotalAccepted1 = 0.0;
    real TotalAccepted2 = 0.0;
    for (int i = 0; i <= Order; i++) {
        if ((op1 == ADD_INTERACTION and op2 == DEL_INTERACTION)
            || (op1 == ADD_DELTA_INTERACTION and op2 == DEL_DELTA_INTERACTION)) {
            if (i == Order)
                continue;
            if (!Equal(Accepted[op1][i] + Accepted[op2][i + 1], 0.0)) {
                TotalAccepted1 += Accepted[op1][i];
                TotalAccepted2 += Accepted[op2][i + 1];
                sprintf(temp, "\t%8s%2i:%15g%15g%15g\n",
                        "Order", i, Accepted[op1][i], Accepted[op2][i + 1],
                        fabs(Accepted[op1][i] - Accepted[op2][i + 1]) / sqrt(Accepted[op1][i] + Accepted[op2][i + 1]));
                Output += temp;
            }
        }
        else if (op1 == JUMP_TO_ORDER0 and op2 == JUMP_BACK_TO_ORDER1) {
            if (i != 0)
                continue;
            if (!Equal(Accepted[op1][i + 1] + Accepted[op2][i], 0.0)) {
                TotalAccepted1 += Accepted[op1][i + 1];
                TotalAccepted2 += Accepted[op2][i];
                sprintf(temp, "\t%8s%2i:%15g%15g%15g\n",
                        "Order", i, Accepted[op1][i + 1], Accepted[op2][i],
                        fabs(Accepted[op1][i + 1] - Accepted[op2][i]) / sqrt(Accepted[op1][i + 1] + Accepted[op2][i]));
                Output += temp;
            }
        }
        else {
            if (!Equal(Accepted[op1][i] + Accepted[op2][i], 0.0)) {
                TotalAccepted1 += Accepted[op1][i];
                TotalAccepted2 += Accepted[op2][i];
                sprintf(temp, "\t%8s%2i:%15g%15g%15g\n",
                        "Order", i, Accepted[op1][i], Accepted[op2][i],
                        fabs(Accepted[op1][i] - Accepted[op2][i]) / sqrt(Accepted[op1][i] + Accepted[op2][i]));
                Output += temp;
            }
        }
    }
    if (!Equal(TotalAccepted1 + TotalAccepted2, 0.0)) {
        sprintf(temp, "\t%10s:%15g%15g%15g\n",
                "Summation", TotalAccepted1, TotalAccepted2,
                fabs(TotalAccepted1 - TotalAccepted2) / sqrt(TotalAccepted1 + TotalAccepted2));
        Output += temp;
    }
    else
        Output += "\tNone\n";
    return Output;
}

void Markov::PrintDetailBalanceInfo()
{
    string Output = "";
    Output = string(60, '=') + "\n";
    Output += "DiagCounter: " + ToString(*Counter) + "\n";
    Output += _DetailBalanceStr(CREATE_WORM);
    Output += _DetailBalanceStr(DELETE_WORM);
    Output += _DetailBalanceStr(MOVE_WORM_G);
    Output += _DetailBalanceStr(MOVE_WORM_W);
    Output += _DetailBalanceStr(RECONNECT);
    Output += _DetailBalanceStr(ADD_INTERACTION);
    Output += _DetailBalanceStr(DEL_INTERACTION);
    Output += _DetailBalanceStr(ADD_DELTA_INTERACTION);
    Output += _DetailBalanceStr(DEL_DELTA_INTERACTION);
    Output += _DetailBalanceStr(CHANGE_TAU_VERTEX);
    Output += _DetailBalanceStr(CHANGE_R_VERTEX);
    Output += _DetailBalanceStr(CHANGE_R_LOOP);
    Output += _DetailBalanceStr(CHANGE_MEASURE_G2W);
    Output += _DetailBalanceStr(CHANGE_MEASURE_W2G);
    Output += _DetailBalanceStr(CHANGE_CONTINUS2DELTA);
    Output += _DetailBalanceStr(CHANGE_DELTA2CONTINUS);
    Output += _DetailBalanceStr(CHANGE_SPIN_VERTEX);
    Output += _DetailBalanceStr(JUMP_TO_ORDER0);
    Output += _DetailBalanceStr(JUMP_BACK_TO_ORDER1);
    Output += string(60, '-') + "\n";
    //    Output += _CheckBalance(CREATE_WORM, DELETE_WORM);
    Output += _CheckBalance(ADD_INTERACTION, DEL_INTERACTION);
    Output += _CheckBalance(ADD_DELTA_INTERACTION, DEL_DELTA_INTERACTION);
    Output += _CheckBalance(CHANGE_MEASURE_G2W, CHANGE_MEASURE_W2G);
    Output += _CheckBalance(CHANGE_CONTINUS2DELTA, CHANGE_DELTA2CONTINUS);
    //    Output += _CheckBalance(JUMP_TO_ORDER0, JUMP_BACK_TO_ORDER1);
    Output += string(60, '=') + "\n";
    LOG_INFO(Output);
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
        if (x < SumofProbofCall[CREATE_WORM])
            CreateWorm();
        //            ;
        else if (x < SumofProbofCall[DELETE_WORM])
            DeleteWorm();
        //            ;
        else if (x < SumofProbofCall[MOVE_WORM_G])
            MoveWormOnG();
        //            ;
        else if (x < SumofProbofCall[MOVE_WORM_W])
            MoveWormOnW();
        //            ;
        else if (x < SumofProbofCall[RECONNECT])
            Reconnect();
        //            ;
        else if (x < SumofProbofCall[ADD_INTERACTION])
            AddInteraction();
        //            ;
        else if (x < SumofProbofCall[DEL_INTERACTION])
            DeleteInteraction();
        //            ;
        else if (x < SumofProbofCall[ADD_DELTA_INTERACTION])
            AddDeltaInteraction();
        //            ;
        else if (x < SumofProbofCall[DEL_DELTA_INTERACTION])
            DeleteDeltaInteraction();
        //            ;
        else if (x < SumofProbofCall[CHANGE_TAU_VERTEX])
            ChangeTauOnVertex();
        //            ;
        else if (x < SumofProbofCall[CHANGE_R_VERTEX])
            ;
        //                    ChangeROnVertex();
        else if (x < SumofProbofCall[CHANGE_R_LOOP])
            ChangeRLoop();
        //            ;
        else if (x < SumofProbofCall[CHANGE_MEASURE_G2W])
            ChangeMeasureFromGToW();
        //            ;
        else if (x < SumofProbofCall[CHANGE_MEASURE_W2G])
            ChangeMeasureFromWToG();
        //            ;
        else if (x < SumofProbofCall[CHANGE_DELTA2CONTINUS])
            ChangeDeltaToContinuous();
        //            ;
        else if (x < SumofProbofCall[CHANGE_CONTINUS2DELTA])
            ChangeContinuousToDelta();
        //            ;
        else if (x < SumofProbofCall[CHANGE_SPIN_VERTEX])
            ;
        //            ChangeSpinOnVertex();
        else if (x < SumofProbofCall[JUMP_TO_ORDER0])
            JumpToOrder0();
        //        ;
        else if (x < SumofProbofCall[JUMP_BACK_TO_ORDER1])
            JumpBackToOrder1();
        //        ;

        (*Counter)++;
    }
}

/**
 *  Create Ira and Masha on a wline
 */
void Markov::CreateWorm()
{
    if (Diag->Order == 0 || Worm->Exist)
        return;

    wLine w = Diag->W.RandomPick(*RNG);
    vertex vin = w->NeighVer(IN);
    vertex vout = w->NeighVer(OUT);

    Momentum kWorm = RandomPickK();
    Momentum kW = w->K - kWorm;
    if (Diag->WHashCheck(kW))
        return;

    int dspin = RandomPickDeltaSpin();
    if (CanNotMoveWorm(dspin, vin->Spin(IN), vin->Spin(OUT)) && CanNotMoveWorm(-dspin, vout->Spin(IN), vout->Spin(OUT)))
        return;

    Complex wWeight = W->Weight(vin->R, vout->R, vin->Tau, vout->Tau,
                                vin->Spin(), vout->Spin(),
                                true, //IsWorm
                                w->IsMeasure, w->IsDelta);

    Complex weightRatio = wWeight / w->Weight;
    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    real wormWeight = weight::Worm::Weight(vin->R, vout->R, vin->Tau, vout->Tau);

    prob *= ProbofCall[DELETE_WORM] / ProbofCall[CREATE_WORM] * (*WormSpaceReweight) * wormWeight * Diag->Order * 2.0;

    Proposed[CREATE_WORM][Diag->Order] += 1.0;
    if (prob >= 1.0 || RNG->urn() < prob) {
        Accepted[CREATE_WORM][Diag->Order] += 1.0;
        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;

        Worm->Exist = true;
        Worm->Ira = vin;
        Worm->Masha = vout;
        Worm->dSpin = dspin;
        Worm->K = kWorm;
        Worm->Weight = wormWeight;

        Diag->ReplaceWHash(w->K, kW);
        w->K = kW;
        w->IsWorm = true;
        w->Weight = wWeight;
    }
}

/**
 *  Delete Ira and Masha on the same wline
 */
void Markov::DeleteWorm()
{
    if (Diag->Order == 0 || !Worm->Exist)
        return;
    vertex &Ira = Worm->Ira;
    vertex &Masha = Worm->Masha;

    wLine w = Ira->NeighW();
    if (!(w == Masha->NeighW()))
        return;
    Momentum k = w->K + SIGN(Ira->Dir) * Worm->K;
    if (Diag->WHashCheck(k))
        return;

    Complex wWeight = W->Weight(Ira->Dir, Ira->R, Masha->R, Ira->Tau, Masha->Tau,
                                Ira->Spin(), Masha->Spin(),
                                false, //IsWorm
                                w->IsMeasure, w->IsDelta);

    Complex weightRatio = wWeight / w->Weight;
    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    prob *= ProbofCall[CREATE_WORM] / (ProbofCall[DELETE_WORM] * (*WormSpaceReweight) * Worm->Weight * Diag->Order * 2.0);

    Proposed[DELETE_WORM][Diag->Order] += 1.0;
    if (prob >= 1.0 || RNG->urn() < prob) {
        Accepted[DELETE_WORM][Diag->Order] += 1.0;
        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;

        Worm->Exist = false;

        w->IsWorm = false;
        Diag->ReplaceWHash(w->K, k);
        w->K = k;
        w->Weight = wWeight;
    }
}

/**
 *  Move Ira along a GLine
 */
void Markov::MoveWormOnG()
{
    if (Diag->Order == 0 || !Worm->Exist)
        return;

    vertex &Ira = Worm->Ira;
    vertex &Masha = Worm->Masha;

    int dir = RandomPickDir();
    gLine g = Ira->NeighG(dir);
    vertex v2 = g->NeighVer(dir);

    if (v2 == Ira || v2 == Masha || CanNotMoveWorm(Worm->dSpin, g->Spin(), dir))
        return;
    Momentum k = g->K - SIGN(dir) * Worm->K;
    if (Diag->GHashCheck(k))
        return;

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
                                 true, //IsWorm
                                 w2->IsMeasure, w2->IsDelta);

    Complex gWeight = G->Weight(INVERSE(dir), Ira->R, v2->R, Ira->Tau, v2->Tau,
                                spinV1[dir], spinV2[INVERSE(dir)], g->IsMeasure);

    Complex weightRatio = w1Weight * w2Weight * gWeight / (g->Weight * w1->Weight * w2->Weight);
    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    real wormWeight = weight::Worm::Weight(v2->R, Masha->R, v2->Tau, Masha->Tau);

    prob *= wormWeight / Worm->Weight;

    Proposed[MOVE_WORM_G][Diag->Order] += 1.0;
    if (prob >= 1.0 || RNG->urn() < prob) {
        Accepted[MOVE_WORM_G][Diag->Order] += 1.0;
        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;

        g->Weight = gWeight;
        Diag->ReplaceGHash(g->K, k);
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
    if (Diag->Order == 0 || !Worm->Exist)
        return;

    vertex &Ira = Worm->Ira;
    vertex &Masha = Worm->Masha;

    wLine w = Ira->NeighW();
    vertex v2 = w->NeighVer(INVERSE(Ira->Dir));
    if (v2 == Ira || v2 == Masha)
        return;
    Momentum k = w->K + SIGN(Ira->Dir) * Worm->K;
    if (Diag->WHashCheck(k))
        return;

    Complex wWeight = W->Weight(Ira->Dir, Ira->R, v2->R, Ira->Tau, v2->Tau, Ira->Spin(),
                                v2->Spin(), w->IsWorm, w->IsMeasure, w->IsDelta);

    Complex weightRatio = wWeight / w->Weight;
    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    real wormWeight = weight::Worm::Weight(v2->R, Masha->R, v2->Tau, Masha->Tau);
    prob *= wormWeight / Worm->Weight;

    Proposed[MOVE_WORM_W][Diag->Order] += 1.0;
    if (prob >= 1.0 || RNG->urn() < prob) {
        Accepted[MOVE_WORM_W][Diag->Order] += 1.0;
        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;

        w->Weight = wWeight;
        Diag->ReplaceWHash(w->K, k);
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
    if (Diag->Order == 0 || !Worm->Exist)
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

    vertex vA = GIA->NeighVer(dir);
    Complex GIAWeight = G->Weight(INVERSE(dir), Masha->R, vA->R, Masha->Tau, vA->Tau,
                                  Masha->Spin(dir), vA->Spin(INVERSE(dir)), GIA->IsMeasure);

    vertex vB = GMB->NeighVer(dir);
    Complex GMBWeight = G->Weight(INVERSE(dir), Ira->R, vB->R, Ira->Tau, vB->Tau,
                                  Ira->Spin(dir), vB->Spin(INVERSE(dir)), GMB->IsMeasure);

    Complex weightRatio = (-1) * GIAWeight * GMBWeight / (GIA->Weight * GMB->Weight);
    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    Proposed[RECONNECT][Diag->Order] += 1.0;
    if (prob >= 1.0 || RNG->urn() < prob) {
        Accepted[RECONNECT][Diag->Order] += 1.0;
        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;
        Diag->SignFermiLoop *= -1;

        Worm->K = k;

        GIA->Weight = GIAWeight;
        GMB->Weight = GMBWeight;

        Masha->nG[dir] = GIA;
        GIA->nVer[INVERSE(dir)] = Masha;

        Ira->nG[dir] = GMB;
        GMB->nVer[INVERSE(dir)] = Ira;
    }
}

/**
 *  add interaction
 *  I---C  ===>  I---A---C
 *  M---D  ===>  M---B---D
 */
void Markov::AddInteraction()
{
    if (!Worm->Exist)
        return;
    if (Diag->Order == 0 || Diag->Order >= Order)
        return;
    vertex Ira = Worm->Ira, Masha = Worm->Masha;

    Momentum kW = RandomPickK();
    if (Diag->WHashCheck(kW))
        return;

    int dir = RandomPickDir();
    int dirW = RandomPickDir();
    gLine GIC = Ira->NeighG(dir), GMD = Masha->NeighG(dir);

    Momentum kIA = GIC->K - SIGN(dir) * SIGN(dirW) * kW;
    Momentum kMB = GMD->K + SIGN(dir) * SIGN(dirW) * kW;
    Momentum kWorm = Worm->K - SIGN(dirW) * kW;
    if (Diag->GHashCheck(kIA))
        return;
    if (Diag->GHashCheck(kMB))
        return;
    if (kIA == kMB)
        return;

    real tauA = RandomPickTau(), tauB = RandomPickTau();

    spin spinA[2] = {GIC->Spin(), GIC->Spin()};
    spin spinB[2] = {GMD->Spin(), GMD->Spin()};
    vertex vC = GIC->NeighVer(dir), vD = GMD->NeighVer(dir);
    Site RA = vC->R, RB = vD->R;

    Complex wWeight = W->Weight(dirW, RA, RB, tauA, tauB, spinA, spinB,
                                false,  //IsWorm
                                false,  //IsMeasure
                                false); //IsDelta

    Complex GIAWeight = G->Weight(INVERSE(dir), Ira->R, RA, Ira->Tau, tauA,
                                  Ira->Spin(dir), spinA[INVERSE(dir)],
                                  false); //IsMeasure

    Complex GMBWeight = G->Weight(INVERSE(dir), Masha->R, RB, Masha->Tau, tauB,
                                  Masha->Spin(dir), spinB[INVERSE(dir)],
                                  false); //IsMeasure

    Complex GACWeight = G->Weight(INVERSE(dir), RA, vC->R, tauA, vC->Tau,
                                  spinA[dir], vC->Spin(INVERSE(dir)), GIC->IsMeasure);

    Complex GBDWeight = G->Weight(INVERSE(dir), RB, vD->R, tauB, vD->Tau,
                                  spinB[dir], vD->Spin(INVERSE(dir)), GMD->IsMeasure);

    Complex weightRatio = (-1) * GIAWeight * GMBWeight * wWeight * GACWeight * GBDWeight / (GIC->Weight * GMD->Weight);

    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    prob *= OrderReWeight[Diag->Order + 1] * ProbofCall[DEL_INTERACTION] / (ProbofCall[ADD_INTERACTION] * OrderReWeight[Diag->Order] * ProbTau(tauA) * ProbTau(tauB));

    Proposed[ADD_INTERACTION][Diag->Order] += 1.0;
    if (prob >= 1.0 || RNG->urn() < prob) {
        Accepted[ADD_INTERACTION][Diag->Order] += 1.0;
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
                      false); //IsMeasure
        Diag->AddGHash(kIA);

        GMB->nVer[INVERSE(dir)] = Masha;
        GMB->nVer[dir] = vB;
        GMB->SetGLine(kMB, GMBWeight,
                      false); //IsMeasure
        Diag->AddGHash(kMB);

        WAB->nVer[dirW] = vA;
        WAB->nVer[INVERSE(dirW)] = vB;
        WAB->SetWLine(kW, wWeight,
                      false,  //IsWorm
                      false,  //IsMeasure
                      false); //IsDelta

        Diag->AddWHash(kW);

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
    if (!Worm->Exist)
        return;
    if (Diag->Order <= 1)
        return;
    vertex Ira = Worm->Ira, Masha = Worm->Masha;

    int dir = RandomPickDir();
    gLine GIA = Ira->NeighG(dir), GMB = Masha->NeighG(dir);
    if (GIA->IsMeasure)
        return;
    if (GMB->IsMeasure)
        return;

    vertex vA = GIA->NeighVer(dir), vB = GMB->NeighVer(dir);
    if (vA->Spin(IN) != vA->Spin(OUT))
        return;
    if (vB->Spin(IN) != vB->Spin(OUT))
        return;

    if (vA->NeighW() != vB->NeighW())
        return;

    wLine wAB = vA->NeighW();
    if (wAB->IsMeasure)
        return;
    if (wAB->IsWorm)
        return;
    if (wAB->IsDelta)
        return;

    gLine GAC = vA->NeighG(dir), GBD = vB->NeighG(dir);
    vertex vC = GAC->NeighVer(dir), vD = GBD->NeighVer(dir);
    if (vA->R != vC->R)
        return;
    if (vB->R != vD->R)
        return;

    Momentum kWorm = Worm->K + SIGN(vA->Dir) * wAB->K;

    Complex GICWeight = G->Weight(INVERSE(dir), Ira->R, vC->R, Ira->Tau, vC->Tau,
                                  Ira->Spin(dir), vC->Spin(INVERSE(dir)), GAC->IsMeasure);

    Complex GMDWeight = G->Weight(INVERSE(dir), Masha->R, vD->R, Masha->Tau, vD->Tau,
                                  Masha->Spin(dir), vD->Spin(INVERSE(dir)), GBD->IsMeasure);

    Complex weightRatio = (-1) * GICWeight * GMDWeight / (GIA->Weight * GMB->Weight * GAC->Weight * GBD->Weight * wAB->Weight);

    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    prob *= OrderReWeight[Diag->Order - 1] * ProbofCall[ADD_INTERACTION] * ProbTau(vA->Tau) * ProbTau(vB->Tau) / (ProbofCall[DEL_INTERACTION] * OrderReWeight[Diag->Order]);

    Proposed[DEL_INTERACTION][Diag->Order] += 1.0;
    if (prob >= 1.0 || RNG->urn() < prob) {
        Accepted[DEL_INTERACTION][Diag->Order] += 1.0;

        Diag->Order--;
        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;

        Diag->Ver.Remove(vA);
        Diag->Ver.Remove(vB);

        Diag->RemoveGHash(GIA->K);
        Diag->G.Remove(GIA);

        Diag->RemoveGHash(GMB->K);
        Diag->G.Remove(GMB);

        Diag->RemoveWHash(wAB->K);
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
 *  add delta interaction
 *  I---C  ===>  I---A---C
 *  M---D  ===>  M---B---D
 */
void Markov::AddDeltaInteraction()
{
    if (!Worm->Exist)
        return;
    if (Diag->Order == 0 || Diag->Order >= Order)
        return;
    vertex Ira = Worm->Ira, Masha = Worm->Masha;

    Momentum kW = RandomPickK();
    if (Diag->WHashCheck(kW))
        return;

    int dir = RandomPickDir();
    int dirW = RandomPickDir();
    gLine GIC = Ira->NeighG(dir), GMD = Masha->NeighG(dir);

    Momentum kIA = GIC->K - SIGN(dir) * SIGN(dirW) * kW;
    Momentum kMB = GMD->K + SIGN(dir) * SIGN(dirW) * kW;
    Momentum kWorm = Worm->K - SIGN(dirW) * kW;
    if (Diag->GHashCheck(kIA))
        return;
    if (Diag->GHashCheck(kMB))
        return;
    if (kIA == kMB)
        return;

    real tauA = RandomPickTau();

    spin spinA[2] = {GIC->Spin(), GIC->Spin()};
    spin spinB[2] = {GMD->Spin(), GMD->Spin()};
    vertex vC = GIC->NeighVer(dir), vD = GMD->NeighVer(dir);
    Site RA = vC->R, RB = vD->R;

    Complex wWeight = W->Weight(dirW, RA, RB, tauA, tauA, spinA, spinB,
                                false, //IsWorm
                                false, //IsMeasure
                                true); //IsDelta

    Complex GIAWeight = G->Weight(INVERSE(dir), Ira->R, RA, Ira->Tau, tauA,
                                  Ira->Spin(dir), spinA[INVERSE(dir)],
                                  false); //IsMeasure

    Complex GMBWeight = G->Weight(INVERSE(dir), Masha->R, RB, Masha->Tau, tauA,
                                  Masha->Spin(dir), spinB[INVERSE(dir)],
                                  false); //IsMeasure

    Complex GACWeight = G->Weight(INVERSE(dir), RA, vC->R, tauA, vC->Tau,
                                  spinA[dir], vC->Spin(INVERSE(dir)), GIC->IsMeasure);

    Complex GBDWeight = G->Weight(INVERSE(dir), RB, vD->R, tauA, vD->Tau,
                                  spinB[dir], vD->Spin(INVERSE(dir)), GMD->IsMeasure);

    Complex weightRatio = (-1) * GIAWeight * GMBWeight * wWeight * GACWeight * GBDWeight / (GIC->Weight * GMD->Weight);

    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    prob *= OrderReWeight[Diag->Order + 1] * ProbofCall[DEL_DELTA_INTERACTION] / (ProbofCall[ADD_DELTA_INTERACTION] * OrderReWeight[Diag->Order] * ProbTau(tauA));

    Proposed[ADD_DELTA_INTERACTION][Diag->Order] += 1.0;
    if (prob >= 1.0 || RNG->urn() < prob) {
        Accepted[ADD_DELTA_INTERACTION][Diag->Order] += 1.0;
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
        vB->SetVertex(RB, tauA, spinB, INVERSE(dirW));

        GIA->nVer[INVERSE(dir)] = Ira;
        GIA->nVer[dir] = vA;
        GIA->SetGLine(kIA, GIAWeight,
                      false); //IsMeasure
        Diag->AddGHash(kIA);

        GMB->nVer[INVERSE(dir)] = Masha;
        GMB->nVer[dir] = vB;
        GMB->SetGLine(kMB, GMBWeight,
                      false); //IsMeasure
        Diag->AddGHash(kMB);

        WAB->nVer[dirW] = vA;
        WAB->nVer[INVERSE(dirW)] = vB;
        WAB->SetWLine(kW, wWeight,
                      false, //IsWorm
                      false, //IsMeasure
                      true); //IsDelta

        Diag->AddWHash(kW);

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
void Markov::DeleteDeltaInteraction()
{
    if (!Worm->Exist)
        return;
    if (Diag->Order <= 1)
        return;
    vertex Ira = Worm->Ira, Masha = Worm->Masha;

    int dir = RandomPickDir();
    gLine GIA = Ira->NeighG(dir), GMB = Masha->NeighG(dir);
    if (GIA->IsMeasure)
        return;
    if (GMB->IsMeasure)
        return;

    vertex vA = GIA->NeighVer(dir), vB = GMB->NeighVer(dir);
    if (vA->Spin(IN) != vA->Spin(OUT))
        return;
    if (vB->Spin(IN) != vB->Spin(OUT))
        return;

    if (vA->NeighW() != vB->NeighW())
        return;

    wLine wAB = vA->NeighW();
    if (wAB->IsMeasure)
        return;
    if (wAB->IsWorm)
        return;
    if (!wAB->IsDelta)
        return;

    gLine GAC = vA->NeighG(dir), GBD = vB->NeighG(dir);
    vertex vC = GAC->NeighVer(dir), vD = GBD->NeighVer(dir);
    if (vA->R != vC->R)
        return;
    if (vB->R != vD->R)
        return;

    Momentum kWorm = Worm->K + SIGN(vA->Dir) * wAB->K;

    Complex GICWeight = G->Weight(INVERSE(dir), Ira->R, vC->R, Ira->Tau, vC->Tau,
                                  Ira->Spin(dir), vC->Spin(INVERSE(dir)), GAC->IsMeasure);

    Complex GMDWeight = G->Weight(INVERSE(dir), Masha->R, vD->R, Masha->Tau, vD->Tau,
                                  Masha->Spin(dir), vD->Spin(INVERSE(dir)), GBD->IsMeasure);

    Complex weightRatio = (-1) * GICWeight * GMDWeight / (GIA->Weight * GMB->Weight * GAC->Weight * GBD->Weight * wAB->Weight);

    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    prob *= OrderReWeight[Diag->Order - 1] * ProbofCall[ADD_DELTA_INTERACTION] * ProbTau(vA->Tau) / (ProbofCall[DEL_DELTA_INTERACTION] * OrderReWeight[Diag->Order]);

    Proposed[DEL_DELTA_INTERACTION][Diag->Order] += 1.0;
    if (prob >= 1.0 || RNG->urn() < prob) {
        Accepted[DEL_DELTA_INTERACTION][Diag->Order] += 1.0;

        Diag->Order--;
        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;

        Diag->Ver.Remove(vA);
        Diag->Ver.Remove(vB);

        Diag->RemoveGHash(GIA->K);
        Diag->G.Remove(GIA);

        Diag->RemoveGHash(GMB->K);
        Diag->G.Remove(GMB);

        Diag->RemoveWHash(wAB->K);
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
    if (Diag->Order == 0 || Worm->Exist)
        return;
    vertex ver = Diag->Ver.RandomPick(*RNG);
    wLine w = ver->NeighW();
    if (w->IsDelta)
        return;

    real tau = RandomPickTau();

    gLine gin = ver->NeighG(IN), gout = ver->NeighG(OUT);
    Complex ginWeight, goutWeight;
    if (gin == gout) {
        //TODO:change to G(-0)
        ginWeight = G->Weight(gin->NeighVer(IN)->R, ver->R,
                              tau, tau,
                              gin->NeighVer(IN)->Spin(OUT), ver->Spin(IN),
                              gin->IsMeasure);
    }
    else {
        ginWeight = G->Weight(gin->NeighVer(IN)->R, ver->R,
                              gin->NeighVer(IN)->Tau, tau,
                              gin->NeighVer(IN)->Spin(OUT), ver->Spin(IN),
                              gin->IsMeasure);
        goutWeight = G->Weight(ver->R, gout->NeighVer(OUT)->R,
                               tau, gout->NeighVer(OUT)->Tau,
                               ver->Spin(OUT), gout->NeighVer(OUT)->Spin(IN),
                               gout->IsMeasure);
    }

    vertex vW = w->NeighVer(INVERSE(ver->Dir));
    Complex wWeight;
    if (vW == ver)
        wWeight = W->Weight(ver->Dir, ver->R, vW->R, tau, tau, ver->Spin(), vW->Spin(),
                            w->IsWorm, w->IsMeasure, w->IsDelta);
    else
        wWeight = W->Weight(ver->Dir, ver->R, vW->R, tau, vW->Tau, ver->Spin(), vW->Spin(),
                            w->IsWorm, w->IsMeasure, w->IsDelta);

    Complex weightRatio = ginWeight * goutWeight * wWeight / (gin->Weight * gout->Weight * w->Weight);

    if (gin == gout)
        weightRatio = ginWeight * wWeight / (gin->Weight * w->Weight);

    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    prob *= ProbTau(ver->Tau) / ProbTau(tau);

    Proposed[CHANGE_TAU_VERTEX][Diag->Order] += 1.0;
    if (prob >= 1.0 || RNG->urn() < prob) {
        Accepted[CHANGE_TAU_VERTEX][Diag->Order] += 1.0;

        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;

        ver->Tau = tau;

        gin->Weight = ginWeight;
        if (gout != gin)
            gout->Weight = goutWeight;
        w->Weight = wWeight;
    }
}

/**
 *  change Spin in one of the two directions of a vertex
 */
void Markov::ChangeSpinOnVertex()
{
    //TODO: If W is spin conserved, return;
    if (Diag->Order == 0 || Worm->Exist)
        return;
    vertex v1 = Diag->Ver.RandomPick(*RNG);
    wLine w1 = v1->NeighW();
    int dir = RandomPickDir();
    gLine g = v1->NeighG(dir);
    vertex v2 = g->NeighVer(dir);
    wLine w2 = v2->NeighW();

    spin spinv1[2] = {v1->Spin(IN), v1->Spin(OUT)};
    spinv1[dir] = FLIP(spinv1[dir]);

    spin spinv2[2] = {v2->Spin(IN), v2->Spin(OUT)};
    spinv2[INVERSE(dir)] = FLIP(spinv2[INVERSE(dir)]);

    spin spinv[2];
    spinv[dir] = spinv1[dir];
    spinv[INVERSE(dir)] = spinv2[INVERSE(dir)];

    Complex gWeight, w1Weight, w2Weight;
    if (w1 == w2) {

        w1Weight = W->Weight(v1->Dir, v1->R, w1->NeighVer(INVERSE(v1->Dir))->R,
                             v1->Tau, w1->NeighVer(INVERSE(v1->Dir))->Tau,
                             spinv, w1->NeighVer(INVERSE(v1->Dir))->Spin(),
                             w1->IsWorm, w1->IsMeasure, w1->IsDelta);
        gWeight = G->Weight(dir, v2->R, v1->R, v2->Tau, v1->Tau,
                            spinv[INVERSE(dir)], spinv[dir], g->IsMeasure);
    }
    else {
        w1Weight = W->Weight(v1->Dir, v1->R, w1->NeighVer(INVERSE(v1->Dir))->R,
                             v1->Tau, w1->NeighVer(INVERSE(v1->Dir))->Tau,
                             spinv1, w1->NeighVer(INVERSE(v1->Dir))->Spin(),
                             w1->IsWorm, w1->IsMeasure, w1->IsDelta);

        w2Weight = W->Weight(v2->Dir, v2->R, w2->NeighVer(INVERSE(v2->Dir))->R,
                             v2->Tau, w2->NeighVer(INVERSE(v2->Dir))->Tau,
                             spinv2, w2->NeighVer(INVERSE(v2->Dir))->Spin(),
                             w2->IsWorm, w2->IsMeasure, w2->IsDelta);

        gWeight = G->Weight(dir, v2->R, v1->R, v2->Tau, v1->Tau,
                            spinv2[INVERSE(dir)], spinv1[dir], g->IsMeasure);
    }

    Complex weightRatio = gWeight * w1Weight * w2Weight / (g->Weight * w1->Weight * w2->Weight);

    if (w1 == w2)
        weightRatio = gWeight * w1Weight / (g->Weight * w1->Weight);

    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    Proposed[CHANGE_SPIN_VERTEX][Diag->Order] += 1.0;
    if (prob >= 1.0 || RNG->urn() < prob) {
        Accepted[CHANGE_SPIN_VERTEX][Diag->Order] += 1.0;
        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;

        v1->SetSpin(spinv1);
        v2->SetSpin(spinv2);
        if (v1 == v2)
            v1->SetSpin(spinv);

        g->Weight = gWeight;
        w1->Weight = w1Weight;
        if (w2 != w1)
            w2->Weight = w2Weight;
    }
}

/**
 *  change R in a vertex
 */
void Markov::ChangeROnVertex()
{
    if (Diag->Order == 0 || Worm->Exist)
        return;
    //TODO: Return if G is local
    vertex ver = Diag->Ver.RandomPick(*RNG);
    Site site = RandomPickSite();
    gLine gin = ver->NeighG(IN), gout = ver->NeighG(OUT);

    Complex ginWeight, goutWeight, wWeight;
    if (gin == gout) {
        ginWeight = G->Weight(site, site,
                              gin->NeighVer(IN)->Tau, ver->Tau,
                              gin->NeighVer(IN)->Spin(OUT), ver->Spin(IN),
                              gin->IsMeasure);
    }
    else {
        ginWeight = G->Weight(gin->NeighVer(IN)->R, site,
                              gin->NeighVer(IN)->Tau, ver->Tau,
                              gin->NeighVer(IN)->Spin(OUT), ver->Spin(IN),
                              gin->IsMeasure);
        goutWeight = G->Weight(site, gout->NeighVer(OUT)->R,
                               ver->Tau, gout->NeighVer(OUT)->Tau,
                               ver->Spin(OUT), gout->NeighVer(OUT)->Spin(IN),
                               gout->IsMeasure);
    }

    wLine w = ver->NeighW();
    vertex vW = w->NeighVer(INVERSE(ver->Dir));

    if (vW == ver)
        wWeight = W->Weight(ver->Dir, site, site, ver->Tau, vW->Tau, ver->Spin(), vW->Spin(),
                            w->IsWorm, w->IsMeasure, w->IsDelta);
    else
        wWeight = W->Weight(ver->Dir, site, vW->R, ver->Tau, vW->Tau, ver->Spin(), vW->Spin(),
                            w->IsWorm, w->IsMeasure, w->IsDelta);

    Complex weightRatio = ginWeight * goutWeight * wWeight / (gin->Weight * gout->Weight * w->Weight);

    if (gin == gout)
        weightRatio = ginWeight * wWeight / (gin->Weight * w->Weight);

    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    prob *= ProbSite(ver->R) / ProbSite(site);

    Proposed[CHANGE_R_VERTEX][Diag->Order] += 1.0;
    if (prob >= 1.0 || RNG->urn() < prob) {
        Accepted[CHANGE_R_VERTEX][Diag->Order] += 1.0;
        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;

        ver->R = site;

        gin->Weight = ginWeight;
        if (gout != gin)
            gout->Weight = goutWeight;
        w->Weight = wWeight;
    }
}

/**
 *  change all R in a fermi loop
 */
void Markov::ChangeRLoop()
{
    if (Diag->Order == 0 || Worm->Exist)
        return;
    //TODO: If G is not a local function, return;
    //TODO: use key word 'static' here to save time
    ASSERT_ALLWAYS(Order <= MAX_ORDER, "Order is too high!");
    vertex v[2 * MAX_ORDER] = {nullptr};
    bool flagVer[2 * MAX_ORDER] = {false};
    int flagW[MAX_ORDER] = {0};
    int n = 0;
    v[0] = Diag->Ver.RandomPick(*RNG);
    Site oldR = v[0]->R;

    while (!flagVer[v[n]->Name]) {
        flagVer[v[n]->Name] = true;
        flagW[v[n]->NeighW()->Name]++;
        v[n + 1] = v[n]->NeighG(OUT)->NeighVer(OUT);

        if (v[n + 1]->R != oldR)
            return;
        n++;
    }

    Site newR = RandomPickSite();

    gLine g = nullptr;
    wLine w = nullptr;

    //TODO: use key word 'static' here to save time
    Complex GWeight[2 * MAX_ORDER] = {Complex(1.0, 0.0)};
    Complex WWeight[2 * MAX_ORDER] = {Complex(1.0, 0.0)};

    Complex oldWeight(1.0, 0.0);
    Complex newWeight(1.0, 0.0);

    for (int i = 0; i < n; i++) {
        g = v[i]->NeighG(OUT);
        GWeight[i] = G->Weight(newR, newR, v[i]->Tau, g->NeighVer(OUT)->Tau,
                               v[i]->Spin(OUT), g->NeighVer(OUT)->Spin(IN), g->IsMeasure);
        newWeight *= GWeight[i];
        oldWeight *= g->Weight;

        w = v[i]->NeighW();
        if (flagW[w->Name] == 1) {
            WWeight[i] = W->Weight(v[i]->Dir, newR, w->NeighVer(INVERSE(v[i]->Dir))->R,
                                   v[i]->Tau, w->NeighVer(INVERSE(v[i]->Dir))->Tau,
                                   v[i]->Spin(), w->NeighVer(INVERSE(v[i]->Dir))->Spin(),
                                   w->IsWorm, w->IsMeasure, w->IsDelta);
            newWeight *= WWeight[i];
            oldWeight *= w->Weight;
        }
        else if (flagW[w->Name] == 2) {
            flagW[w->Name] = 0;
            WWeight[i] = W->Weight(v[i]->Dir, newR, newR,
                                   v[i]->Tau, w->NeighVer(INVERSE(v[i]->Dir))->Tau,
                                   v[i]->Spin(), w->NeighVer(INVERSE(v[i]->Dir))->Spin(),
                                   w->IsWorm, w->IsMeasure, w->IsDelta);
            newWeight *= WWeight[i];
            oldWeight *= w->Weight;
        }
        else if (flagW[w->Name] == 0) {
            WWeight[i] = W->Weight(v[i]->Dir, newR, newR,
                                   v[i]->Tau, w->NeighVer(INVERSE(v[i]->Dir))->Tau,
                                   v[i]->Spin(), w->NeighVer(INVERSE(v[i]->Dir))->Spin(),
                                   w->IsWorm, w->IsMeasure, w->IsDelta);
        }
    }

    Complex weightRatio = newWeight / oldWeight;
    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    prob *= ProbSite(oldR) / ProbSite(newR);

    Proposed[CHANGE_R_LOOP][Diag->Order] += 1.0;
    if (prob >= 1.0 || RNG->urn() < prob) {
        Accepted[CHANGE_R_LOOP][Diag->Order] += 1.0;
        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;

        for (int i = 0; i < n; i++) {
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
    if (Diag->Order == 0 || Worm->Exist || !Diag->MeasureGLine)
        return;

    wLine w = Diag->W.RandomPick(*RNG);
    if (w->IsDelta)
        return;

    gLine g = Diag->GMeasure;
    Complex gWeight = G->Weight(g->NeighVer(IN)->R, g->NeighVer(OUT)->R,
                                g->NeighVer(IN)->Tau, g->NeighVer(OUT)->Tau,
                                g->Spin(), g->Spin(),
                                false); //IsMeasure

    Complex wWeight = W->Weight(w->NeighVer(IN)->R, w->NeighVer(OUT)->R,
                                w->NeighVer(IN)->Tau, w->NeighVer(OUT)->Tau,
                                w->NeighVer(IN)->Spin(), w->NeighVer(OUT)->Spin(),
                                w->IsWorm,
                                true, //IsMeasure
                                w->IsDelta);

    Complex weightRatio = gWeight * wWeight / (g->Weight * w->Weight);
    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    //proposal probility: (1/2N)/(1/N)
    prob *= 0.5 * (*PolarReweight) * ProbofCall[CHANGE_MEASURE_W2G] / (ProbofCall[CHANGE_MEASURE_G2W]);

    Proposed[CHANGE_MEASURE_G2W][Diag->Order] += 1.0;
    if (prob >= 1.0 || RNG->urn() < prob) {
        Accepted[CHANGE_MEASURE_G2W][Diag->Order] += 1.0;
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
    if (Diag->Order == 0 || Worm->Exist || Diag->MeasureGLine)
        return;

    gLine g = Diag->G.RandomPick(*RNG);

    wLine w = Diag->WMeasure;
    if (w->IsDelta)
        return;

    Complex gWeight = G->Weight(g->NeighVer(IN)->R, g->NeighVer(OUT)->R,
                                g->NeighVer(IN)->Tau, g->NeighVer(OUT)->Tau,
                                g->Spin(), g->Spin(),
                                true); //IsMeasure

    Complex wWeight = W->Weight(w->NeighVer(IN)->R, w->NeighVer(OUT)->R,
                                w->NeighVer(IN)->Tau, w->NeighVer(OUT)->Tau,
                                w->NeighVer(IN)->Spin(), w->NeighVer(OUT)->Spin(),
                                w->IsWorm,
                                false, //IsMeasure
                                w->IsDelta);

    Complex weightRatio = gWeight * wWeight / (g->Weight * w->Weight);
    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    prob *= ProbofCall[CHANGE_MEASURE_G2W] / (0.5 * (*PolarReweight) * ProbofCall[CHANGE_MEASURE_W2G]);

    Proposed[CHANGE_MEASURE_W2G][Diag->Order] += 1.0;
    if (prob >= 1.0 || RNG->urn() < prob) {
        Accepted[CHANGE_MEASURE_W2G][Diag->Order] += 1.0;
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
    if (Diag->Order < 2 || Worm->Exist)
        return;
    wLine w = Diag->W.RandomPick(*RNG);
    if ((!w->IsDelta) || w->IsMeasure)
        return;
    vertex vin = w->NeighVer(IN), vout = w->NeighVer(OUT);
    gLine G1 = vout->NeighG(IN), G2 = vout->NeighG(OUT);
    real tau = RandomPickTau();
    Complex wWeight = W->Weight(vin->R, vout->R, vin->Tau, tau, vin->Spin(),
                                vout->Spin(), w->IsWorm, w->IsMeasure,
                                false); //IsDelta

    Complex G1Weight, G2Weight, weightRatio;
    if (G1 == G2) {
        G1Weight = G->Weight(G1->NeighVer(IN)->R, vout->R,
                             tau, tau,
                             G1->NeighVer(IN)->Spin(OUT), vout->Spin(IN),
                             G1->IsMeasure);
        weightRatio = G1Weight * wWeight / (G1->Weight * w->Weight);
    }
    else {
        G1Weight = G->Weight(G1->NeighVer(IN)->R, vout->R,
                             G1->NeighVer(IN)->Tau, tau,
                             G1->NeighVer(IN)->Spin(OUT), vout->Spin(IN),
                             G1->IsMeasure);

        G2Weight = G->Weight(OUT, G2->NeighVer(OUT)->R, vout->R,
                             G2->NeighVer(OUT)->Tau, tau,
                             G2->NeighVer(OUT)->Spin(IN), vout->Spin(OUT),
                             G2->IsMeasure);
        weightRatio = G1Weight * G2Weight * wWeight / (G1->Weight * G2->Weight * w->Weight);
    }

    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    prob *= ProbofCall[CHANGE_CONTINUS2DELTA] / (ProbofCall[CHANGE_DELTA2CONTINUS] * ProbTau(tau));

    Proposed[CHANGE_DELTA2CONTINUS][Diag->Order] += 1.0;
    if (prob >= 1.0 || RNG->urn() < prob) {
        Accepted[CHANGE_DELTA2CONTINUS][Diag->Order] += 1.0;
        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;

        vout->Tau = tau;
        w->IsDelta = false;
        w->Weight = wWeight;

        G1->Weight = G1Weight;
        if (G1 != G2)
            G2->Weight = G2Weight;
    }
}

/**
 *  Change a Wline from continuous to Delta
 */
void Markov::ChangeContinuousToDelta()
{
    if (Diag->Order < 2 || Worm->Exist)
        return;

    wLine w = Diag->W.RandomPick(*RNG);
    if (w->IsDelta || w->IsMeasure)
        return;

    vertex vin = w->NeighVer(IN), vout = w->NeighVer(OUT);
    gLine G1 = vout->NeighG(IN), G2 = vout->NeighG(OUT);

    Complex wWeight = W->Weight(vin->R, vout->R, vin->Tau, vin->Tau, vin->Spin(),
                                vout->Spin(), w->IsWorm, w->IsMeasure,
                                true); //IsDelta

    Complex G1Weight, G2Weight, weightRatio;
    if (G1 == G2) {
        G1Weight = G->Weight(G1->NeighVer(IN)->R, vout->R,
                             vin->Tau, vin->Tau,
                             G1->NeighVer(IN)->Spin(OUT), vout->Spin(IN),
                             G1->IsMeasure);
        weightRatio = G1Weight * wWeight / (G1->Weight * w->Weight);
    }
    else {
        G1Weight = G->Weight(G1->NeighVer(IN)->R, vout->R,
                             G1->NeighVer(IN)->Tau, vin->Tau,
                             G1->NeighVer(IN)->Spin(OUT), vout->Spin(IN),
                             G1->IsMeasure);
        G2Weight = G->Weight(OUT, G2->NeighVer(OUT)->R, vout->R,
                             G2->NeighVer(OUT)->Tau, vin->Tau,
                             G2->NeighVer(OUT)->Spin(IN), vout->Spin(OUT),
                             G2->IsMeasure);
        weightRatio = G1Weight * G2Weight * wWeight / (G1->Weight * G2->Weight * w->Weight);
    }

    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    prob *= ProbofCall[CHANGE_DELTA2CONTINUS] * ProbTau(vout->Tau) / ProbofCall[CHANGE_CONTINUS2DELTA];

    Proposed[CHANGE_CONTINUS2DELTA][Diag->Order] += 1.0;
    if (prob >= 1.0 || RNG->urn() < prob) {
        Accepted[CHANGE_CONTINUS2DELTA][Diag->Order] += 1.0;
        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;

        vout->Tau = vin->Tau;
        w->IsDelta = true;
        w->Weight = wWeight;

        G1->Weight = G1Weight;
        if (G1 != G2)
            G2->Weight = G2Weight;
    }
}

void Markov::JumpToOrder0()
{
    if (Worm->Exist || Diag->Order != 1)
        return;

    vertex Ver1 = &Diag->Ver[0];
    vertex Ver2 = &Diag->Ver[1];
    if (Ver1->R != Ver2->R)
        return;

    Complex weightRatio;
    weightRatio = weight::Norm::Weight() / Diag->Weight;

    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    prob *= (ProbofCall[JUMP_BACK_TO_ORDER1] * ProbSite(Ver1->R) * ProbTau(Ver1->Tau) * ProbTau(Ver2->Tau) * 0.5 * 0.5 * OrderReWeight[0]) / (ProbofCall[JUMP_TO_ORDER0] * OrderReWeight[1]);

    Proposed[JUMP_TO_ORDER0][Diag->Order] += 1.0;
    if (prob >= 1.0 || RNG->urn() < prob) {
        Accepted[JUMP_TO_ORDER0][Diag->Order] += 1.0;
        Diag->Order = 0;
        Diag->Phase *= sgn;
        Diag->Weight = weight::Norm::Weight();
    }
}

void Markov::JumpBackToOrder1()
{
    if (Worm->Exist || Diag->Order != 0)
        return;

    Site R = RandomPickSite();
    real Tau1 = RandomPickTau();
    real Tau2 = RandomPickTau();
    spin SpinG1 = RandomPickSpin();
    spin SpinG2 = RandomPickSpin();

    vertex Ver1 = &Diag->Ver[0];
    vertex Ver2 = &Diag->Ver[1];
    gLine G1 = Ver1->NeighG(OUT);
    gLine G2 = Ver2->NeighG(OUT);
    wLine W1 = Ver1->NeighW();

    spin SpinV1[2] = {SpinG2, SpinG1};
    spin SpinV2[2] = {SpinG1, SpinG2};

    Complex weightG1 = G->Weight(R, R, Tau1, Tau2, SpinG1, SpinG1, G1->IsMeasure);
    Complex weightG2 = G->Weight(R, R, Tau2, Tau1, SpinG2, SpinG2, G2->IsMeasure);
    Complex weightW = W->Weight(Ver1->Dir, R, R, Tau1, Tau2, SpinV1, SpinV2, false, W1->IsMeasure, W1->IsDelta);

    Complex weightRatio = -1.0 * weightG1 * weightG2 * weightW / Diag->Weight;
    real prob = mod(weightRatio);
    Complex sgn = phase(weightRatio);

    prob *= ProbofCall[JUMP_TO_ORDER0] * OrderReWeight[1] / (ProbofCall[JUMP_BACK_TO_ORDER1] * OrderReWeight[0] * ProbSite(R) * ProbTau(Tau1) * ProbTau(Tau2) * 0.5 * 0.5);

    Proposed[JUMP_BACK_TO_ORDER1][Diag->Order] += 1.0;
    if (prob >= 1.0 || RNG->urn() < prob) {
        Accepted[JUMP_BACK_TO_ORDER1][Diag->Order] += 1.0;
        Diag->Order = 1;

        Diag->Phase *= sgn;
        Diag->Weight *= weightRatio;

        Ver1->R = R;
        Ver2->R = R;
        Ver1->Tau = Tau1;
        Ver2->Tau = Tau2;
        Ver1->SetSpin(SpinV1);
        Ver2->SetSpin(SpinV2);

        G1->Weight = weightG1;
        G2->Weight = weightG2;
        W1->Weight = weightW;
    }
}

Momentum Markov::RandomPickK()
{
    return (Momentum)(RNG->irn(-MAX_K, MAX_K));
}

int Markov::RandomPickDeltaSpin()
{
    return RNG->irn(0, 1) * 2 - 1;
}

spin Markov::RandomPickSpin()
{
    return (RNG->irn(0, 1) == 0 ? DOWN : UP);
}

int Markov::RandomPickDir()
{
    return RNG->irn(0, 1);
}

real Markov::RandomPickTau()
{
    return RNG->urn() * Beta;
}

real Markov::ProbTau(real tau)
{
    return 1.0 / Beta;
}

bool Markov::RandomPickBool()
{
    return (RNG->irn(0, 1) == 0 ? true : false);
}

Site Markov::RandomPickSite()
{
    Vec<int> coord;
    for (int i = 0; i < D; i++)
        coord[i] = RNG->irn(0, Lat->Size[i] - 1);
    return (Site(RNG->irn(0, Lat->SublatVol - 1), coord));
}

real Markov::ProbSite(const Site &site)
{

    return 1.0 / (Lat->Vol * Lat->SublatVol);
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
