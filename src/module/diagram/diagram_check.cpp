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
bool Diagram::CheckTopo()
{
    if(G.HowMany()!=2*Order)
        ABORT("Number of G is wrong!");
    if(W.HowMany()!=Order)
        ABORT("Number of W is wrong!");
    if(Ver.HowMany()!=2*Order)
        ABORT("Number of Vertex is wrong!");
    
    for (int i = 0; i < G.HowMany(); i++) {
        for (int dir = 0; i < 2; i++) {
            vertex v = G(i)->NeighVer(dir);
            if (!Ver.Exist(v))
                ABORT("nVer not exists!" + v->PrettyString());
            if (G(i) != v->NeighG(INVERSE(dir)))
                ABORT("Neigh of G is incorrect!" + G(i)->PrettyString());
        }
    }
    
    for (int i = 0; i < W.HowMany(); i++) {
        for (int dir = 0; i < 2; i++) {
            vertex v = W(i)->NeighVer(dir);
            if (!Ver.Exist(v))
                ABORT("nVer not exists!" + v->PrettyString());
            if (W(i) != v->NeighW())
                ABORT("Neigh of W is incorrect!" + W(i)->PrettyString());
            if (v->Dir!=dir)
                ABORT("Direction of Vertex is incorrect!" + v->PrettyString());
        }
    }
    
    return true;
}

bool Diagram::CheckStatus()
{
    int totalmeasure=0;
    for( int i=0; i< G.HowMany(); i++)
    {
        if(G(i)->IsMeasure)
            totalmeasure +=1;
    }
    if(totalmeasure!=(MeasureGLine?1:0))
        ABORT("number of Measuring Gline is wrong!");
    
    totalmeasure = 0;
    for( int i=0; i< W.HowMany(); i++)
    {
        if(W(i)->IsMeasure)
            totalmeasure +=1;
    }
    if(totalmeasure!=(MeasureGLine?0:1))
        ABORT("number of Measuring Wline is wrong!");
    
    return true;
}

bool Diagram::CheckK()
{
    return true;
}

bool Diagram::CheckSpin()
{
    return true;
}

bool Diagram::CheckSite()
{
    return true;
}

bool Diagram::CheckTau()
{
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
    if (!CheckTopo())
        return false;
    if (!CheckWeight())
        return false;
    return true;
}
