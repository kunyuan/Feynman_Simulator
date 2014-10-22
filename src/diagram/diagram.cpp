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
#include "weight.h"
#include "rng.h"
using namespace std;
bool Diagram::IsWorm(const Vertex &v)
{
    if (Worm.Ira == v.Name || Worm.Masha == v.Name)
        return true;
    else
        return false;
}

const string ToString(const spin &s)
{
    if (s == UP)
        return "UP";
    else
        return "DOWN";
}

Diagram::Diagram()
    : Order(0), Phase(Complex(1.0, 0.0)), Weight(Complex(1.0, 0.0)), G("GLine"), W("WLine"), Ver("Vertex")
{
}

/****************   GLine  *****************************/
spin Diagram::Spin(GLine &g, const int &dir)
{
    return Ver[g.Vertex[dir]].Spin[FlipDir(dir)];
}

int Diagram::Sublattice(GLine &g, const int &dir)
{
    return Ver[g.Vertex[dir]].R.Sublattice;
}

string Diagram::PrettyString(GLine &g)
{
    stringstream os;
    os << "{V " << g.Vertex[IN] << "}->-" << ToString(Spin(g, IN)) << "---";
    os << "[G " << g.Name << " ,K:" << g.K << ",Weight:" << g.Weight << "]";
    os << "---" << ToString(Spin(g, OUT)) << "->-{V " << g.Vertex[OUT] << "}";
    os << endl;
    return os.str();
}

/****************   WLine  *****************************/
spin Diagram::Spin(WLine &w, const int &dir1, const int &dir2)
{
    return Ver[w.Vertex[dir1]].Spin[dir2];
}

int Diagram::Sublattice(WLine &w, const int &dir)
{
    return Ver[w.Vertex[dir]].R.Sublattice;
}

string Diagram::PrettyString(WLine &w)
{
    stringstream os;
    os << "{V " << w.Vertex[IN] << "| " << ToString(Spin(w, IN, IN)) << "," << ToString(Spin(w, IN, OUT)) << "}~~~";
    os << "<W " << w.Name << ", K:" << w.K << ",Weight:" << w.Weight << ">";
    os << "~~~"
       << "{" << ToString(Spin(w, OUT, IN)) << "," << ToString(Spin(w, OUT, OUT)) << "|W " << w.Vertex[OUT] << "}";
    os << endl;
    return os.str();
}

/****************   Vertex  *****************************/
spin Diagram::Spin(Vertex &v, const int &dir)
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
    os << "-[G " << v.G[IN] << "]-" << ToString(v.Spin[IN]) << "->--";
    os << "{V " << v.Name << ",r:" << v.R.Coordinate.PrettyString() << ",tau:" << v.Tau << "}";
    os << "-->-" << ToString(v.Spin[OUT]) << "-[G " << v.G[OUT] << "]-";
    os << "  /~~~<W " << v.W << ">";
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

GLine &Diagram::RandomPickG()
{
    return G[RNG.irn(0, G.HowMany())];
}

WLine &Diagram::RandomPickW()
{
    return W[RNG.irn(0, W.HowMany())];
}

Vertex &Diagram::RandomPickVer()
{
    return Ver[RNG.irn(0, Ver.HowMany())];
}

bool Diagram::FixDiagram()
{
    Order = W.HowMany();
    for (int index = 0; index < G.HowMany(); index++) {
        GLine &g = G[index];
        for (int dir = 0; dir < 2; dir++) {
            NeighVer(g, dir).G[FlipDir(dir)] = index;
        }
    }

    for (int index = 0; index < W.HowMany(); index++) {
        WLine &w = W[index];
        w.IsWorm = false;

        for (int dir = 0; dir < 2; dir++) {
            NeighVer(w, dir).W = index;
        }
    }

    for (int index = 0; index < Ver.HowMany(); index++) {
        //TODO: Do something here if you want to fix vertex
    }
    
    Phase = Complex(1.0, 0.0);
    Weight = Complex(1.0, 0.0);
    return true;
}
