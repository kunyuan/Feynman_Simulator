//
//  diagram_global_check.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/10/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "diagram.h"
#include "utility/abort.h"
#include "module/weight/weight.h"

using namespace std;
using namespace diag;

///*************************   Diagram check    *************************/
bool Diagram::CheckG()
{
    for (int i = 0; i < G.HowMany(); i++) {
        for (int dir = 0; i < 2; i++) {
            vertex v = G(i)->NeighVer(dir);
            if (!Ver.Exist(v))
                ABORT("nVer not exists!" + v->PrettyString());
            if (G(i) != v->NeighG(INVERSE(dir)))
                ABORT("Neigh of G is incorrect!" + G(i)->PrettyString());
        }
    }
    return true;
}

bool Diagram::CheckW()
{
    for (int i = 0; i < W.HowMany(); i++) {
        for (int dir = 0; i < 2; i++) {
            vertex v = W(i)->NeighVer(dir);
            if (!Ver.Exist(v))
                ABORT("nVer not exists!" + v->PrettyString());
            if (W(i) != v->NeighW())
                ABORT("Neigh of W is incorrect!" + W(i)->PrettyString());
        }
    }
    return true;
}

bool Diagram::CheckVer()
{
    //TODO: check more with vertex if you want!
    return true;
}

bool Diagram::CheckWeight()
{
    if (DEBUGMODE && (GWeight == nullptr || WWeight == nullptr))
        ABORT("G and W weight are not defined yet!");

    Complex DiagWeight(1.0, 0.0);
    Complex gWeight, wWeight;
    vertex vin, vout;

    for (int i = 0; i < G.HowMany(); i++) {
        DiagWeight *= G(i)->Weight;

        vin = G(i)->NeighVer(IN);
        vout = G(i)->NeighVer(OUT);
        gWeight = GWeight->Weight(vin->R, vout->R, vin->Tau, vout->Tau,
                                  vin->Spin(OUT), vout->Spin(IN), G(i)->IsMeasure);
        if (!Equal(G(i)->Weight, gWeight))
            return false;
    }
    for (int i = 0; i < W.HowMany(); i++) {
        DiagWeight *= W(i)->Weight;

        vin = W(i)->NeighVer(IN);
        vout = W(i)->NeighVer(OUT);
        wWeight = WWeight->Weight(vin->R, vout->R, vin->Tau, vout->Tau, vin->Spin(),
                                  vout->Spin(), W(i)->IsWorm, W(i)->IsMeasure, W(i)->IsDelta);
        if (!Equal(W(i)->Weight, wWeight))
            return false;
    }
    DiagWeight *= SignFermiLoop * (Order % 2 == 0 ? 1 : -1);

    if (!Equal(DiagWeight, Weight))
        return false;
    return true;
}

bool Diagram::CheckDiagram()
{
    //TODO: don't forget to check diagram weight
    if (!CheckG())
        return false;
    if (!CheckW())
        return false;
    if (!CheckVer())
        return false;
    if (!CheckWeight())
        return false;
    return true;
}
