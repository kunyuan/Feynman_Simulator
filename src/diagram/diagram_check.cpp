//
//  diagram_global_check.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/10/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "diagram.h"
#include "../utility/abort.h"
using namespace std;

///*************************   Diagram check    *************************/
bool Diagram::CheckG()
{
    for (int i = 0; i < G.HowMany(); i++) {
        for (int dir = 0; i < 2; i++) {
            Vertex &v = NeighVer(G[i], dir);
            if (!Ver.Exist(v))
                ABORT("Vertex not exists!" + PrettyString(v));
            if (&G[i] != &NeighG(v, 1 - dir))
                ABORT("Neigh of G is incorrect!" + PrettyString(G[i]));
        }
    }
    return true;
}

bool Diagram::CheckW()
{
    for (int i = 0; i < W.HowMany(); i++) {
        for (int dir = 0; i < 2; i++) {
            Vertex &v = NeighVer(W[i], dir);
            if (!Ver.Exist(v))
                ABORT("Vertex not exists!" + PrettyString(v));
            if (&W[i] != &NeighW(v))
                ABORT("Neigh of W is incorrect!" + PrettyString(W[i]));
        }
    }
    return true;
}

bool Diagram::CheckVer()
{
    //TODO: check more with vertex if you want!
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
    return true;
}
