//
//  diagram_global.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/9/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "diagram.h"
#include "../observable/weight.h"
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
    Lat = nullptr;
    GWeight = nullptr;
    WWeight = nullptr;
}

void Diagram::Build(Lattice *lat, Weight::G *g, Weight::W *w)
{
    Lat = lat;
    GWeight = g;
    WWeight = w;
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
    return G[RNG.irn(0, G.HowMany() - 1)];
}

WLine &Diagram::RandomPickW()
{
    return W[RNG.irn(0, W.HowMany() - 1)];
}

Vertex &Diagram::RandomPickVer()
{
    return Ver[RNG.irn(0, Ver.HowMany() - 1)];
}

void Diagram::ClearDiagram()
{
    while (G.HowMany() > 0)
        G.Remove(G.HowMany() - 1);
    while (W.HowMany() > 0)
        W.Remove(W.HowMany() - 1);
    while (Ver.HowMany() > 0)
        Ver.Remove(Ver.HowMany() - 1);
}
bool Diagram::FixDiagram()
{
    if (DEBUGMODE && Lat == nullptr)
        ABORT("Lattice is not defined yet!");
    if (DEBUGMODE && (GWeight == nullptr || WWeight == nullptr))
        ABORT("G and W weight are not defined yet!");

    Order = W.HowMany();
    Worm.Exist = false;

    Weight = Complex(1.0, 0.0);
    for (int index = 0; index < G.HowMany(); index++) {
        GLine &g = G[index];
        Vertex &vin = NeighVer(g, IN);
        Vertex &vout = NeighVer(g, OUT);

        vin.G[OUT] = index;
        vout.G[IN] = index;

        g.Weight = GWeight->Weight(Lat->Distance(vin.R, vout.R), vout.Tau - vin.Tau, vin.Spin[OUT], vout.Spin[IN]);
        Weight *= g.Weight;
    }

    for (int index = 0; index < W.HowMany(); index++) {
        WLine &w = W[index];
        Vertex &vin = NeighVer(w, IN);
        Vertex &vout = NeighVer(w, OUT);

        w.IsWorm = false;

        vin.W = index;
        vout.W = index;

        w.Weight = WWeight->Weight(Lat->Distance(vin.R, vout.R), vout.Tau - vin.Tau, vin.Spin, vout.Spin, w.IsWorm);
        Weight *= w.Weight;
    }

    for (int index = 0; index < Ver.HowMany(); index++) {
        //TODO: Do something here if you want to fix vertex
    }

    Phase = phase(Weight);
    return true;
}
