//
//  diagram_global_check.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/10/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "diagram.h"
#include "abort.h"
using namespace std;

///*************************   Diagram check    *************************/
bool Diagram::CheckG()
{
    for (int i = 0; i < G.HowMany(); i++) {
        for (int dir = 0; i < 2; i++) {
            Vertex &v = NeighVer(G[i], dir);
            if (!Ver.Exist(v))
                ABORT("Vertex not exists!" + v.PrettyString());
            if (&G[i] != &NeighG(v, 1 - dir))
                ABORT("Neigh of G is incorrect!" + G[i].PrettyString());
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
                ABORT("Vertex not exists!" + v.PrettyString());
            if (&W[i] != &NeighW(v))
                ABORT("Neigh of W is incorrect!" + W[i].PrettyString());
            if (W[i].Spin[dir][IN] != VerSpin(v, IN))
                ABORT("Spin of W in direction " + to_string(dir) + ",0 is incorrect!" + W[i].PrettyString());
            if (W[i].Spin[dir][OUT] != VerSpin(v, OUT))
                ABORT("Spin of W in direction " + to_string(dir) + ",1 is incorrect!" + W[i].PrettyString());
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
