//
// Created by yuan huang on 10/5/17.
//
#include "diag_calculator.h"
#include "utility/dictionary.h"
#include <iostream>
#include <stdio.h>
#include "math.h"
#include "utility/utility.h"
#include "utility/momentum.h"
#include "module/diagram/diagram.h"
#include "module/parameter/parameter.h"
#include "lattice/lattice.h"
#include "module/weight/weight.h"
#include "module/weight/component.h"
#include <array>

using namespace diagCalc;
using namespace std;
using namespace para;

#define NAME(x) #x

int Find_Position(std::array<int, 2*MAX_ORDER> arr, int value);

bool MonteCarlo::BuildNew(ParaMC &para, DiagramDict &diag, weight::Weight &weight)
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

    OperationName[CHANGE_R] = NAME(CHANGE_R);
    OperationName[CHANGE_TAU] = NAME(CHANGE_TAU);
    OperationName[CHANGE_DELTA2CONTINUS] = NAME(CHANGE_DELTA2CONTINUS);
    OperationName[CHANGE_CONTINUS2DELTA] = NAME(CHANGE_CONTINUS2DELTA);
    OperationName[CHANGE_MEASURE_G2W] = NAME(CHANGE_MEASURE_G2W);
    OperationName[CHANGE_MEASURE_W2G] = NAME(CHANGE_MEASURE_W2G);
    OperationName[ADD_INTERACTION] = NAME(ADD_INTERACTION);
    OperationName[DEL_INTERACTION] = NAME(DEL_INTERACTION);
    OperationName[JUMP_TO_ORDER0] = NAME(JUMP_TO_ORDER0);
    OperationName[JUMP_TO_ORDER1] = NAME(JUMP_TO_ORDER1);

    New_DiagConf.Initialize();

    IndexToSpin = {UP, DOWN};

    InitializeGWTable();
    Old_G_Weight = New_G_Weight;
    Old_W_Weight = New_W_Weight;

    Complex initial_weight = SumAllDiagrams();

    New_DiagConf.Weight = mod(initial_weight);
    New_DiagConf.Phase = phase(initial_weight);

    Old_DiagConf = New_DiagConf;

    return true;
}


void MonteCarlo::Reset(ParaMC &para, DiagramDict &diag, weight::Weight &weight)
{
    Beta = para.Beta;
    Order = para.Order;
    Counter = &para.Counter;
    Lat = &para.Lat;
    OrderReWeight = para.OrderReWeight.data();
    PolarReweight = &para.PolarReweight;
    Diag = &diag;

    //WormSpaceReweight = &para.WormSpaceReweight;
    //Worm = &diag.Worm;
    Sigma = weight.Sigma;
    Polar = weight.Polar;
    G = weight.G;
    W = weight.W;
    RNG = &para.RNG;
}

void MonteCarlo::InitializeGWTable() {

    for (int in = 0; in < 2*New_DiagConf.Order; in ++){
        for (int out = 0; out < 2*New_DiagConf.Order; out ++) {

            Site RIn = New_DiagConf.R_list[in];
            Site ROut = New_DiagConf.R_list[out];
            real TauIn = New_DiagConf.Tau_list[in];
            real TauOut = New_DiagConf.Tau_list[out];

            for (int sp = 0; sp < 2; sp ++){
                Complex gWeight = G->Weight(IN, RIn, ROut, TauIn, TauOut,
                                            IndexToSpin[sp], IndexToSpin[sp], false);
                New_G_Weight[sp][in][out] = gWeight;
            }

            for(int sp_leftin = 0; sp_leftin < 2; sp_leftin ++){
                for(int sp_leftout = 0; sp_leftout < 2; sp_leftout ++) {
                    for (int sp_rightin = 0; sp_rightin < 2; sp_rightin++) {
                        for (int sp_rightout = 0; sp_rightout < 2; sp_rightout++) {
                            spin SpLeft[2] = {IndexToSpin[sp_leftin], IndexToSpin[sp_leftout]};
                            spin SpRight[2] = {IndexToSpin[sp_rightin], IndexToSpin[sp_rightout]};
                            Complex wWeight = W->Weight(IN, RIn, ROut, TauIn, TauOut,
                                                        SpLeft, SpRight, false, false, false);
                            New_W_Weight[SpinIndex(SpLeft, SpRight)][in][out] = wWeight;
                        }
                    }
                }
            }
        }
    }
}

void MonteCarlo::UpdateGWTable(int vertex) {
    int in = vertex;
    for (int out = 0; out < 2*New_DiagConf.Order; out ++){

        Site RIn = New_DiagConf.R_list[in];
        Site ROut = New_DiagConf.R_list[out];
        real TauIn = New_DiagConf.Tau_list[in];
        real TauOut = New_DiagConf.Tau_list[out];

        for (int sp = 0; sp < 2; sp ++){
            Complex gWeight = G->Weight(IN, RIn, ROut, TauIn, TauOut,
                                        IndexToSpin[sp], IndexToSpin[sp], false);
            New_G_Weight[sp][in][out] = gWeight;
        }

        for(int sp_leftin = 0; sp_leftin < 2; sp_leftin ++){
            for(int sp_leftout = 0; sp_leftout < 2; sp_leftout ++) {
                for (int sp_rightin = 0; sp_rightin < 2; sp_rightin++) {
                    for (int sp_rightout = 0; sp_rightout < 2; sp_rightout++) {
                        spin SpLeft[2] = {IndexToSpin[sp_leftin], IndexToSpin[sp_leftout]};
                        spin SpRight[2] = {IndexToSpin[sp_rightin], IndexToSpin[sp_rightout]};
                        Complex wWeight = W->Weight(IN, RIn, ROut, TauIn, TauOut,
                                                    SpLeft, SpRight, false, false, false);
                        New_W_Weight[SpinIndex(SpLeft, SpRight)][in][out] = wWeight;
                    }
                }
            }
        }
    }

    int out = vertex;
    for (int in = 0; in < 2*New_DiagConf.Order; in ++){

        Site RIn = New_DiagConf.R_list[in];
        Site ROut = New_DiagConf.R_list[out];
        real TauIn = New_DiagConf.Tau_list[in];
        real TauOut = New_DiagConf.Tau_list[out];

        for (int sp = 0; sp < 2; sp ++){
            Complex gWeight = G->Weight(IN, RIn, ROut, TauIn, TauOut,
                                        IndexToSpin[sp], IndexToSpin[sp], false);
            New_G_Weight[sp][in][out] = gWeight;
        }

        for(int sp_leftin = 0; sp_leftin < 2; sp_leftin ++){
            for(int sp_leftout = 0; sp_leftout < 2; sp_leftout ++) {
                for (int sp_rightin = 0; sp_rightin < 2; sp_rightin++) {
                    for (int sp_rightout = 0; sp_rightout < 2; sp_rightout++) {
                        spin SpLeft[2] = {IndexToSpin[sp_leftin], IndexToSpin[sp_leftout]};
                        spin SpRight[2] = {IndexToSpin[sp_rightin], IndexToSpin[sp_rightout]};
                        Complex wWeight = W->Weight(IN, RIn, ROut, TauIn, TauOut,
                                                    SpLeft, SpRight, false, false, false);
                        New_W_Weight[SpinIndex(SpLeft, SpRight)][in][out] = wWeight;
                    }
                }
            }
        }
    }
}

int Find_Position(std::array<int, 2*MAX_ORDER> arr, int value){
    auto iter = std::find(arr.begin(), arr.end(), value);
    return int(std::distance(arr.begin(), iter));
}

Complex MonteCarlo::SumAllDiagrams() {

    Complex total_weight=Complex(0.0, 0.0);

    int order = New_DiagConf.Order;

    DiagramsList permutations = Diag->AllDiagramConfig[order-1];
    SpinsList spin_configs = Diag->AllSpinConfig[order-1];
    FermiSignList fermi_signs = Diag->AllFermiSignConfig[order-1];

    Complex weight;

    for(int i = 0; i < Diag->AllDiagramConfig.size(); i++){
        auto permutation = permutations[i];
        auto spins = spin_configs[i];
        weight = Complex(1.0, 0.0);
        for (int in=0; in < 2*order; in++){
            weight *= New_G_Weight[spins[in]][in][permutation[in]];
        }
        for (int in=0; in < order; in++){

            int leftout =in*2, rightout = in*2+1;
            int leftin = Find_Position(permutation, in*2), rightin = Find_Position(permutation, in*2+1);
            spin SpLeft[2] = {IndexToSpin[spins[leftin]], IndexToSpin[spins[leftout]]};
            spin SpRight[2] = {IndexToSpin[spins[rightin]], IndexToSpin[spins[rightout]]};
            int spin_index = SpinIndex(SpLeft, SpRight);
            weight *= New_W_Weight[spin_index][in*2][in*2+1];
        }

        total_weight += weight * fermi_signs[i];
    }
    return total_weight;
}

std::string MonteCarlo::_DetailBalanceStr(Operations op)
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

std::string MonteCarlo::_CheckBalance(Operations op1, Operations op2)
{
    string Output = OperationName[op1] + "<------>" + OperationName[op2] + ":\n";
    char temp[80];
    real TotalAccepted1 = 0.0;
    real TotalAccepted2 = 0.0;
    for (int i = 0; i <= Order; i++) {
        if (op1 == ADD_INTERACTION and op2 == DEL_INTERACTION){
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
        else if (op1 == JUMP_TO_ORDER0 and op2 == JUMP_TO_ORDER1) {
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

void MonteCarlo::PrintDetailBalanceInfo()
{
    string Output = "";
    Output = string(60, '=') + "\n";
    Output += "DiagCounter: " + ToString(*Counter) + "\n";
    Output += _DetailBalanceStr(CHANGE_R);
    Output += _DetailBalanceStr(CHANGE_TAU);
    Output += _DetailBalanceStr(CHANGE_DELTA2CONTINUS);
    Output += _DetailBalanceStr(CHANGE_CONTINUS2DELTA);
    Output += _DetailBalanceStr(CHANGE_MEASURE_G2W);
    Output += _DetailBalanceStr(CHANGE_MEASURE_W2G);
    Output += _DetailBalanceStr(ADD_INTERACTION);
    Output += _DetailBalanceStr(DEL_INTERACTION);
    Output += _DetailBalanceStr(JUMP_TO_ORDER0);
    Output += _DetailBalanceStr(JUMP_TO_ORDER1);
    Output += string(60, '-') + "\n";
    Output += _CheckBalance(CHANGE_DELTA2CONTINUS, CHANGE_CONTINUS2DELTA);
    Output += _CheckBalance(CHANGE_MEASURE_G2W, CHANGE_MEASURE_W2G);
    Output += _CheckBalance(ADD_INTERACTION, DEL_INTERACTION);
    Output += _CheckBalance(JUMP_TO_ORDER0, JUMP_TO_ORDER1);
    Output += string(60, '=') + "\n";
    LOG_INFO(Output);
}

/**
*  \brief let the Grasshopper hops for Steps
*
*  @param Steps
*/
void MonteCarlo::Hop(int sweep)
{
    for (int i = 0; i < sweep; i++) {
        double x = RNG->urn();
        if (x < SumofProbofCall[CHANGE_R])
            ChangeR();
        else if (x < SumofProbofCall[CHANGE_TAU])
            ChangeTau();
        else if (x < SumofProbofCall[CHANGE_DELTA2CONTINUS])
            ChangeDeltaToContinuous();
        else if (x < SumofProbofCall[CHANGE_CONTINUS2DELTA])
            ChangeContinuousToDelta();
        else if (x < SumofProbofCall[CHANGE_MEASURE_G2W])
            ChangeMeasureFromGToW();
        else if (x < SumofProbofCall[CHANGE_MEASURE_W2G])
            ChangeMeasureFromWToG();
        else if (x < SumofProbofCall[ADD_INTERACTION])
            AddInteraction();
        else if (x < SumofProbofCall[DEL_INTERACTION])
            DeleteInteraction();
        else if (x < SumofProbofCall[JUMP_TO_ORDER0])
            JumpToOrder0();
        else if (x < SumofProbofCall[JUMP_TO_ORDER1])
            JumpToOrder1();

        (*Counter)++;
    }
}


/**
 *  change space variable R
 */
void MonteCarlo::ChangeR() {
    if (Old_DiagConf.Order == 0)
        return;

    int vertex = RandomPickVertex(Old_DiagConf.Order);
    Site site = RandomPickSite();

    New_DiagConf.R_list[vertex] = site;

    UpdateGWTable(vertex);

    Complex new_weight = SumAllDiagrams();

    real prob = mod(new_weight/Old_DiagConf.Weight);
    prob *= ProbSite(Old_DiagConf.R_list[vertex]) / ProbSite(site);

    Proposed[CHANGE_R][Old_DiagConf.Order] += 1.0;
    if (prob >= 1.0 || RNG->urn() < prob) {
        LOG_INFO("change R");
        Accepted[CHANGE_R][Old_DiagConf.Order] += 1.0;

        New_DiagConf.Weight = mod(new_weight);
        New_DiagConf.Phase = phase(new_weight);

        Old_DiagConf = New_DiagConf;
    }
}

/**
 *  change time variable tau
 */
void MonteCarlo::ChangeTau() {
    if (Old_DiagConf.Order == 0)
        return;

    int vertex = RandomPickVertex(Old_DiagConf.Order);
    real tau = RandomPickTau();

    New_DiagConf.Tau_list[vertex] = tau;

    UpdateGWTable(vertex);

    Complex new_weight = SumAllDiagrams();

    real prob = mod(new_weight/Old_DiagConf.Weight);
    prob *= ProbTau(Old_DiagConf.Tau_list[vertex]) / ProbTau(tau);

    Proposed[CHANGE_TAU][Old_DiagConf.Order] += 1.0;
    if (prob >= 1.0 || RNG->urn() < prob) {
        LOG_INFO("change tau");
        Accepted[CHANGE_TAU][Old_DiagConf.Order] += 1.0;

        New_DiagConf.Weight = mod(new_weight);
        New_DiagConf.Phase = phase(new_weight);

        Old_DiagConf = New_DiagConf;
    }
}

/**
 *  change measuring line from G to W
 */
void MonteCarlo::ChangeMeasureFromGToW()
{
//    LOG_INFO("change measure line from g to w");
//    if (Diag->Order == 0 || Worm->Exist || !Diag->MeasureGLine)
//        return;
//
//    wLine w = Diag->W.RandomPick(*RNG);
//    if (w->IsDelta)
//        return;
//
//    gLine g = Diag->GMeasure;
//    Complex gWeight = G->Weight(g->NeighVer(IN)->R, g->NeighVer(OUT)->R,
//                                g->NeighVer(IN)->Tau, g->NeighVer(OUT)->Tau,
//                                g->Spin(), g->Spin(),
//                                false); //IsMeasure
//
//    Complex wWeight = W->Weight(w->NeighVer(IN)->R, w->NeighVer(OUT)->R,
//                                w->NeighVer(IN)->Tau, w->NeighVer(OUT)->Tau,
//                                w->NeighVer(IN)->Spin(), w->NeighVer(OUT)->Spin(),
//                                w->IsWorm,
//                                true, //IsMeasure
//                                w->IsDelta);
//
//    Complex weightRatio = gWeight * wWeight / (g->Weight * w->Weight);
//    real prob = mod(weightRatio);
//    Complex sgn = phase(weightRatio);
//
//    //proposal probility: (1/2N)/(1/N)
//    prob *= 0.5 * (*PolarReweight) * ProbofCall[CHANGE_MEASURE_W2G] / (ProbofCall[CHANGE_MEASURE_G2W]);
//
//    Proposed[CHANGE_MEASURE_G2W][Diag->Order] += 1.0;
//    if (prob >= 1.0 || RNG->urn() < prob) {
//        Accepted[CHANGE_MEASURE_G2W][Diag->Order] += 1.0;
//        Diag->Phase *= sgn;
//        Diag->Weight *= weightRatio;
//
//        Diag->MeasureGLine = false;
//        Diag->GMeasure = nullptr;
//        Diag->WMeasure = w;
//
//        g->IsMeasure = false;
//        w->IsMeasure = true;
//
//        g->Weight = gWeight;
//        w->Weight = wWeight;
//    }
}

/**
 *  change measuring line from W to G
 */
void MonteCarlo::ChangeMeasureFromWToG()
{
//    LOG_INFO("change measure line from w to g");
//    if (Diag->Order == 0 || Worm->Exist || Diag->MeasureGLine)
//        return;
//
//    gLine g = Diag->G.RandomPick(*RNG);
//
//    wLine w = Diag->WMeasure;
//    if (w->IsDelta)
//        return;
//
//    Complex gWeight = G->Weight(g->NeighVer(IN)->R, g->NeighVer(OUT)->R,
//                                g->NeighVer(IN)->Tau, g->NeighVer(OUT)->Tau,
//                                g->Spin(), g->Spin(),
//                                true); //IsMeasure
//
//    Complex wWeight = W->Weight(w->NeighVer(IN)->R, w->NeighVer(OUT)->R,
//                                w->NeighVer(IN)->Tau, w->NeighVer(OUT)->Tau,
//                                w->NeighVer(IN)->Spin(), w->NeighVer(OUT)->Spin(),
//                                w->IsWorm,
//                                false, //IsMeasure
//                                w->IsDelta);
//
//    Complex weightRatio = gWeight * wWeight / (g->Weight * w->Weight);
//    real prob = mod(weightRatio);
//    Complex sgn = phase(weightRatio);
//
//    prob *= ProbofCall[CHANGE_MEASURE_G2W] / (0.5 * (*PolarReweight) * ProbofCall[CHANGE_MEASURE_W2G]);
//
//    Proposed[CHANGE_MEASURE_W2G][Diag->Order] += 1.0;
//    if (prob >= 1.0 || RNG->urn() < prob) {
//        Accepted[CHANGE_MEASURE_W2G][Diag->Order] += 1.0;
//        Diag->Phase *= sgn;
//        Diag->Weight *= weightRatio;
//
//        Diag->MeasureGLine = true;
//        Diag->GMeasure = g;
//        Diag->WMeasure = nullptr;
//
//        g->IsMeasure = true;
//        w->IsMeasure = false;
//
//        g->Weight = gWeight;
//        w->Weight = wWeight;
//    }
}

/**
 *  change a Wline from Delta to continuous
 */
void MonteCarlo::ChangeDeltaToContinuous()
{
//    LOG_INFO("change w line from delta to continuous");
//    if (Diag->Order < 2 || Worm->Exist)
//        return;
//    wLine w = Diag->W.RandomPick(*RNG);
//    if ((!w->IsDelta) || w->IsMeasure)
//        return;
//    vertex vin = w->NeighVer(IN), vout = w->NeighVer(OUT);
//    gLine G1 = vout->NeighG(IN), G2 = vout->NeighG(OUT);
//    real tau = RandomPickTau();
//    Complex wWeight = W->Weight(vin->R, vout->R, vin->Tau, tau, vin->Spin(),
//                                vout->Spin(), w->IsWorm, w->IsMeasure,
//                                false); //IsDelta
//
//    Complex G1Weight, G2Weight, weightRatio;
//    if (G1 == G2) {
//        G1Weight = G->Weight(G1->NeighVer(IN)->R, vout->R,
//                             tau, tau,
//                             G1->NeighVer(IN)->Spin(OUT), vout->Spin(IN),
//                             G1->IsMeasure);
//        weightRatio = G1Weight * wWeight / (G1->Weight * w->Weight);
//    }
//    else {
//        G1Weight = G->Weight(G1->NeighVer(IN)->R, vout->R,
//                             G1->NeighVer(IN)->Tau, tau,
//                             G1->NeighVer(IN)->Spin(OUT), vout->Spin(IN),
//                             G1->IsMeasure);
//
//        G2Weight = G->Weight(OUT, G2->NeighVer(OUT)->R, vout->R,
//                             G2->NeighVer(OUT)->Tau, tau,
//                             G2->NeighVer(OUT)->Spin(IN), vout->Spin(OUT),
//                             G2->IsMeasure);
//        weightRatio = G1Weight * G2Weight * wWeight / (G1->Weight * G2->Weight * w->Weight);
//    }
//
//    real prob = mod(weightRatio);
//    Complex sgn = phase(weightRatio);
//
//    prob *= ProbofCall[CHANGE_CONTINUS2DELTA] / (ProbofCall[CHANGE_DELTA2CONTINUS] * ProbTau(tau));
//
//    Proposed[CHANGE_DELTA2CONTINUS][Diag->Order] += 1.0;
//    if (prob >= 1.0 || RNG->urn() < prob) {
//        Accepted[CHANGE_DELTA2CONTINUS][Diag->Order] += 1.0;
//        Diag->Phase *= sgn;
//        Diag->Weight *= weightRatio;
//
//        vout->Tau = tau;
//        w->IsDelta = false;
//        w->Weight = wWeight;
//
//        G1->Weight = G1Weight;
//        if (G1 != G2)
//            G2->Weight = G2Weight;
//    }
}

/**
 *  Change a Wline from continuous to Delta
 */
void MonteCarlo::ChangeContinuousToDelta()
{
//    LOG_INFO("change w line from continuous to delta");
//    if (Diag->Order < 2 || Worm->Exist)
//        return;
//
//    wLine w = Diag->W.RandomPick(*RNG);
//    if (w->IsDelta || w->IsMeasure)
//        return;
//
//    vertex vin = w->NeighVer(IN), vout = w->NeighVer(OUT);
//    gLine G1 = vout->NeighG(IN), G2 = vout->NeighG(OUT);
//
//    Complex wWeight = W->Weight(vin->R, vout->R, vin->Tau, vin->Tau, vin->Spin(),
//                                vout->Spin(), w->IsWorm, w->IsMeasure,
//                                true); //IsDelta
//
//    Complex G1Weight, G2Weight, weightRatio;
//    if (G1 == G2) {
//        G1Weight = G->Weight(G1->NeighVer(IN)->R, vout->R,
//                             vin->Tau, vin->Tau,
//                             G1->NeighVer(IN)->Spin(OUT), vout->Spin(IN),
//                             G1->IsMeasure);
//        weightRatio = G1Weight * wWeight / (G1->Weight * w->Weight);
//    }
//    else {
//        G1Weight = G->Weight(G1->NeighVer(IN)->R, vout->R,
//                             G1->NeighVer(IN)->Tau, vin->Tau,
//                             G1->NeighVer(IN)->Spin(OUT), vout->Spin(IN),
//                             G1->IsMeasure);
//        G2Weight = G->Weight(OUT, G2->NeighVer(OUT)->R, vout->R,
//                             G2->NeighVer(OUT)->Tau, vin->Tau,
//                             G2->NeighVer(OUT)->Spin(IN), vout->Spin(OUT),
//                             G2->IsMeasure);
//        weightRatio = G1Weight * G2Weight * wWeight / (G1->Weight * G2->Weight * w->Weight);
//    }
//
//    real prob = mod(weightRatio);
//    Complex sgn = phase(weightRatio);
//
//    prob *= ProbofCall[CHANGE_DELTA2CONTINUS] * ProbTau(vout->Tau) / ProbofCall[CHANGE_CONTINUS2DELTA];
//
//    Proposed[CHANGE_CONTINUS2DELTA][Diag->Order] += 1.0;
//    if (prob >= 1.0 || RNG->urn() < prob) {
//        Accepted[CHANGE_CONTINUS2DELTA][Diag->Order] += 1.0;
//        Diag->Phase *= sgn;
//        Diag->Weight *= weightRatio;
//
//        vout->Tau = vin->Tau;
//        w->IsDelta = true;
//        w->Weight = wWeight;
//
//        G1->Weight = G1Weight;
//        if (G1 != G2)
//            G2->Weight = G2Weight;
//    }
}

/**
 *  increase diagram order
 */
void MonteCarlo::AddInteraction() {
//    LOG_INFO("increase diagram order");
}

/**
 *  decrease diagram order
 */
void MonteCarlo::DeleteInteraction() {
//    LOG_INFO("decrease diagram order");
}


/**
 *  jump to order 0 diagram from order 1
 */
void MonteCarlo::JumpToOrder0() {
//    LOG_INFO("jump to order 0 diagram");
}

/**
 *  jump to order 1 diagram from order 0
 */
void MonteCarlo::JumpToOrder1() {
//    LOG_INFO("jump to order 1 diagram");
}

//First In/Out: direction of WLine; Second In/Out: direction of Vertex
int MonteCarlo::SpinIndex(spin SpinInIn, spin SpinInOut, spin SpinOutIn, spin SpinOutOut)
{
    return SpinInIn * SPIN3 + SpinInOut * SPIN2 + SpinOutIn * SPIN + SpinOutOut;
}

int MonteCarlo::SpinIndex(const spin* TwoSpinIn, const spin* TwoSpinOut)
{
    return SpinIndex(TwoSpinIn[0], TwoSpinIn[1],
                     TwoSpinOut[0], TwoSpinOut[1]);
}

int MonteCarlo::RandomPickVertex(int order)
{
    return RNG->irn(0, 2*order);
}

real MonteCarlo::RandomPickTau()
{
    return RNG->urn() * Beta;
}

real MonteCarlo::ProbTau(real tau)
{
    return 1.0 / Beta;
}

Site MonteCarlo::RandomPickSite()
{
    Vec<int> coord;
    for (int i = 0; i < D; i++)
        coord[i] = RNG->irn(0, Lat->Size[i] - 1);
    return (Site(RNG->irn(0, Lat->SublatVol - 1), coord));
}

real MonteCarlo::ProbSite(const Site &site)
{
    return 1.0 / (Lat->Vol * Lat->SublatVol);
}

