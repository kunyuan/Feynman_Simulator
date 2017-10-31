//
//  diagram_global_check.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/10/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "diagram.h"
#include "utility/abort.h"
#include "module/weight/component.h"
#include <iostream>

using namespace std;
using namespace diag;

bool Diagram::CheckDiagram()
{
    if (!DEBUGMODE)
        return true;

    if (MeasureGammaGW!=0)
        return true;

    if (!_CheckTopo())
        return false;
    if (!_CheckStatus())
        return false;
    if (!_CheckK())
        return false;
    if (!_CheckSpin())
        return false;
    if (!_CheckWeight())
        return false;
    return true;
}
///*************************   Diagram check    *************************/
bool Diagram::_CheckTopo()
{
    if(MeasureGammaGW!=0)
        return true;

    if (Order == 0)
        return true;
    if (G.HowMany() != 2 * Order){
        WriteDiagram2gv("error_diagram.gv");
        ABORT("Number of G is wrong!");
    }
    if (W.HowMany() != Order)
        ABORT("Number of W is wrong!");
    if (Ver.HowMany() != 2 * Order)
        ABORT("Number of Vertex is wrong!");

    for (int i = 0; i < G.HowMany(); i++) {
        for (int dir = 0; dir < 2; dir++) {
            vertex v = G(i)->NeighVer(dir);
            if (!Ver.Exist(v))
                ABORT("nVer not exists!" + v->PrettyString());
            if (G(i) != v->NeighG(INVERSE(dir)))
                ABORT("Neigh of G is incorrect!" + G(i)->PrettyString());
        }
    }

    for (int i = 0; i < W.HowMany(); i++) {
        for (int dir = 0; dir < 2; dir++) {
            vertex v = W(i)->NeighVer(dir);
            if (!Ver.Exist(v))
                ABORT("nVer not exists!" + v->PrettyString());
            if (W(i) != v->NeighW())
                ABORT("Neigh of W is incorrect!" + W(i)->PrettyString());
            if (v->Dir != dir)
                ABORT("Direction of Vertex is incorrect!" + v->PrettyString());
        }
    }

    return true;
}

bool Diagram::_CheckStatus()
{
    if (Order == 0)
        return true;
    int totalmeasure = 0;
    for (int i = 0; i < G.HowMany(); i++) {
        if (G(i)->IsMeasure) {
            totalmeasure += 1;
            if (G(i) != GMeasure)
                ABORT("GMeasure error!" + G(i)->PrettyString());
        }
    }
    if (totalmeasure != (MeasureGLine ? 1 : 0))
        ABORT("number of Measuring Gline is wrong!");

    totalmeasure = 0;
    for (int i = 0; i < W.HowMany(); i++) {
        if (W(i)->IsMeasure)
            totalmeasure += 1;
        if (W(i)->IsWorm && Worm.Ira->NeighW() != W(i) && Worm.Masha->NeighW() != W(i))
            ABORT("W IsWorm status error! no worm!" + W(i)->PrettyString());
        if (Worm.Exist)
            if ((!W(i)->IsWorm) && (Worm.Ira->NeighW() == W(i) || Worm.Masha->NeighW() == W(i)))
                ABORT("W IsWorm status error! has worm!" + W(i)->PrettyString());
        if (W(i)->IsDelta) {
            if (W(i)->NeighVer(IN)->Tau != W(i)->NeighVer(OUT)->Tau)
                ABORT("W is delta function, tau error!" + W(i)->PrettyString());
            if (W(i)->IsMeasure)
                ABORT("W is delta function and measuring line, error!" + W(i)->PrettyString());
        }
    }
    if (totalmeasure != (MeasureGLine ? 0 : 1))
        ABORT("number of Measuring Wline is wrong!");

    return true;
}

bool Diagram::_CheckK()
{
    if(MeasureGammaGW!=0)
        return true;
    if (Order == 0)
        return true;
    Momentum totalk;
    for (int i = 0; i < Ver.HowMany(); i++) {
        totalk = Ver(i)->NeighG(IN)->K - Ver(i)->NeighG(OUT)->K;
        if (Ver(i)->Dir == IN)
            totalk -= Ver(i)->NeighW()->K;
        else
            totalk += Ver(i)->NeighW()->K;
        if (Worm.Exist && Ver(i) == Worm.Ira)
            totalk -= Worm.K;
        else if (Worm.Exist && Ver(i) == Worm.Masha)
            totalk += Worm.K;
        if (totalk != 0)
            ABORT("K is not conserved!" + Ver(i)->PrettyString());
    }
    return true;
}

bool Diagram::_CheckSpin()
{
    if (Order == 0)
        return true;

    for (int i = 0; i < G.HowMany(); i++) {
        if (G(i)->NeighVer(IN)->Spin(OUT) != G(i)->NeighVer(OUT)->Spin(IN))
            ABORT("The spin on Gline is not the same" + G(i)->PrettyString());
    }
    //TODO: check Spin on W only when spin is conserved in W
    return true;
}

bool Diagram::_CheckWeight()
{
    cout << "Check weight..." << endl;
    if (Order == 0)
        return Equal(Weight, weight::Norm::Weight());
    else {
        Complex DiagWeight(1.0, 0.0);
        Complex gWeight, wWeight;
        vertex vin, vout;

        for (int i = 0; i < G.HowMany(); i++) {
            DiagWeight *= G(i)->Weight;

            vin = G(i)->NeighVer(IN);
            vout = G(i)->NeighVer(OUT);
            gWeight = GWeight->Weight(vin->R, vout->R, vin->Tau, vout->Tau,
                                      vin->Spin(OUT), vout->Spin(IN), G(i)->IsMeasure, G(i)->IsGammaG, &UExt);
            if (!Equal(G(i)->Weight, gWeight))
                return false;
            if (Equal(G(i)->Weight, Complex(0.0, 0.0)))
                return false;
        }
        for (int i = 0; i < W.HowMany(); i++) {
            DiagWeight *= W(i)->Weight;

            vin = W(i)->NeighVer(IN);
            vout = W(i)->NeighVer(OUT);
            wWeight = WWeight->Weight(vin->R, vout->R, vin->Tau, vout->Tau, vin->Spin(),
                                      vout->Spin(), W(i)->IsWorm, W(i)->IsMeasure, W(i)->IsDelta, W(i)->IsGammaW, &UExt);
            if (!Equal(W(i)->Weight, wWeight))
                return false;
            if (Equal(W(i)->Weight, Complex(0.0, 0.0)))
                return false;
        }
        DiagWeight *= SignFermiLoop * (Order % 2 == 0 ? 1 : -1);
        return Equal(DiagWeight, Weight);
    }
}
