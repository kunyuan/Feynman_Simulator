//
//  diagram_global.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/9/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "diagram.h"
#include "utility.h"
using namespace std;

Diagram::Diagram()
    : G("GLine"), W("WLine"), Ver("Vertex")
{
}

/****************   GLine  *****************************/
spin Diagram::G_Spin(GLine &g, int dir)
{
    return g.Spin[dir];
}

int Diagram::G_Sublattice(GLine &g, int dir)
{
    return Ver[g.Vertex[dir]].Sublattice;
}

/****************   Vertex  *****************************/
spin Diagram::VerSpin(Vertex &v, int dir)
{
    return G[v.G[dir]].Spin[1 - dir];
}

/****************   Diagram  *****************************/

Vertex &Diagram::NeighVer(GLine &g, int dir)
{
    return Ver[g.Vertex[dir]];
}

Vertex &Diagram::NeighVer(WLine &w, int dir)
{
    return Ver[w.Vertex[dir]];
}

GLine &Diagram::NeighG(Vertex &v, int dir)
{
    return G[v.G[dir]];
}

WLine &Diagram::NeighW(Vertex &v)
{
    return W[v.W];
}

bool Diagram::FixDiagram()
{
    //TODO: you may also need to fix diagram weight
    for (int index = 0; index < Ver.HowMany(); index++) {
        for (int dir = 0; dir < 2; dir++) {
            NeighVer(G[index], dir).G[1 - dir] = index;
        }
    }

    for (int index = 0; index < Ver.HowMany(); index++) {
        WLine &w = W[index];
        for (int dir = 0; dir < 2; dir++) {
            Vertex &v = NeighVer(w, dir);
            v.W = index;
            w.Spin[dir][0] = VerSpin(v, 0);
            w.Spin[dir][1] = VerSpin(v, 1);
        }
    }

    for (int index = 0; index < Ver.HowMany(); index++) {
        //TODO: Do something here if you want to fix vertex
    }
    return true;
}
