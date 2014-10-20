//
//  diagram_global.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/9/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "diagram.h"
#include "utility.h"
#include "convention.h"
#include "lattice.h"
using namespace std;

ostream &operator<<(ostream &os, spin &s)
{
    if (s == UP)
        return os << "UP";
    else
        return os << "DOWN";
}

Diagram::Diagram()
    : G("GLine"), W("WLine"), Ver("Vertex")
{
}

/****************   GLine  *****************************/
spin Diagram::Spin(GLine &g, const int& dir)
{
    return Ver[g.Vertex[dir]].Spin[FlipDir(dir)];
}

int Diagram::Sublattice(GLine &g, const int& dir)
{
    return Ver[g.Vertex[dir]].R.Sublattice;
}

string Diagram::PrettyString(GLine &g)
{
    stringstream os;
    os << "\n";
    os << g.Vertex[IN] << " (" << Spin(g,IN) << ") >>===";
    os << "Name:" << g.Name << ",Weight:" << g.Weight;
    os << "===>>" << g.Vertex[OUT] << " (" << Spin(g,OUT) << ");";
    os << endl;
    return os.str();
}



/****************   WLine  *****************************/
spin Diagram::Spin(WLine &w, const int& dir1, const int& dir2)
{
    return Ver[w.Vertex[dir1]].Spin[dir2];
}

int Diagram::Sublattice(WLine &w, const int& dir)
{
    return Ver[w.Vertex[dir]].R.Sublattice;
}


string Diagram::PrettyString(WLine& w)
{
    stringstream os;
    os << "\n";
    os << w.Vertex[IN] << " (" << Spin(w, IN, IN) << "," << Spin(w, IN,OUT) << ") ~~~";
    os << "Name:" << w.Name << ",Weight:" << w.Weight;
    os << "~~~" << w.Vertex[OUT] << " (" << Spin(w,OUT,IN) << "," << Spin(w,OUT,OUT) << ");";
    os << endl;
    return os.str();
}

/****************   Vertex  *****************************/
spin Diagram::Spin(Vertex &v, const int& dir)
{
    return v.Spin[dir];
}

int Diagram::Sublattice(Vertex &v)
{
    return v.R.Sublattice;
}

string Diagram::PrettyString(Vertex &v)
{
    stringstream os;
    os << "\n";
    os << v.G[IN] << ">===";
    os << "Name:" << v.Name << ",r:" << v.R.Coordinate.PrettyString() << ",tau:" << v.tau;
    os << "===>" << v.G[OUT] << "\n";
    os << "       ~~~~" << v.W;
    os << endl;
    return os.str();
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
            NeighVer(G[index], dir).G[FlipDir(dir)] = index;
        }
    }

    for (int index = 0; index < Ver.HowMany(); index++) {
        WLine &w = W[index];
        for (int dir = 0; dir < 2; dir++) {
            Vertex &v = NeighVer(w, dir);
            v.W = index;
        }
    }

    for (int index = 0; index < Ver.HowMany(); index++) {
        //TODO: Do something here if you want to fix vertex
    }
    return true;
}