//
//  diagram_global.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/9/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "diagram.h"
#include "../observable/weight.h"
#include "../utility/rng.h"
using namespace std;
bool Diagram::IsWorm(vertex v)
{
    if (Worm.Ira == v || Worm.Masha == v)
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
    : Order(0), Phase(Complex(1.0, 0.0)), Weight(Complex(1.0, 0.0)), G("GLine"), W("WLine"), Ver("nVer")
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
spin Diagram::Spin(gLine g, int dir)
{
    return g->nVer[dir]->Spin[FlipDir(dir)];
}

int Diagram::Sublattice(gLine g, int dir)
{
    return g->nVer[dir]->R.Sublattice;
}

string Diagram::PrettyString(gLine g)
{
    stringstream os;
    os << "{V " << g->nVer[IN] << "}->-" << ToString(Spin(g, IN)) << "---";
    os << "[G " << g->Name << " ,K:" << g->K << ",Weight:" << g->Weight << "]";
    os << "---" << ToString(Spin(g, OUT)) << "->-{V " << g->nVer[OUT] << "}";
    return os.str();
}

/****************   WLine  *****************************/
spin Diagram::Spin(wLine w, int dir1, int dir2)
{
    return w->nVer[dir1]->Spin[dir2];
}

int Diagram::Sublattice(wLine w, int dir)
{
    return w->nVer[dir]->R.Sublattice;
}

string Diagram::PrettyString(wLine w)
{
    stringstream os;
    os << "{V " << w->nVer[IN] << "| " << ToString(Spin(w, IN, IN)) << "," << ToString(Spin(w, IN, OUT)) << "}~~~";
    os << "<W " << w->Name << ", K:" << w->K << ",Weight:" << w->Weight << ">";
    os << "~~~"
       << "{" << ToString(Spin(w, OUT, IN)) << "," << ToString(Spin(w, OUT, OUT)) << "|W " << w->nVer[OUT] << "}";
    return os.str();
}

/****************   Vertex  *****************************/
spin Diagram::Spin(vertex v, int dir)
{
    return v->Spin[dir];
}

int Diagram::Sublattice(vertex v)
{
    return v->R.Sublattice;
}

string Diagram::PrettyString(vertex v)
{
    stringstream os;
    os << "-[G " << v->nG[IN] << "]-" << ToString(v->Spin[IN]) << "->--";
    os << "{V " << v->Name << ",r:" << v->R.Coordinate.PrettyString() << ",tau:" << v->Tau << "}";
    os << "-->-" << ToString(v->Spin[OUT]) << "-[G " << v->nG[OUT] << "]-";
    os << "  /~~~<W " << v->nW << ">";
    return os.str();
}

/****************   Diagram  *****************************/

vertex Diagram::NeighVer(gLine g, int dir)
{
    return g->nVer[dir];
}

vertex Diagram::NeighVer(wLine w, int dir)
{
    return w->nVer[dir];
}

gLine Diagram::NeighG(vertex v, int dir)
{
    return v->nG[dir];
}

wLine Diagram::NeighW(vertex v)
{
    return v->nW;
}

gLine Diagram::RandomPickG()
{
    return &G[RNG.irn(0, G.HowMany() - 1)];
}

wLine Diagram::RandomPickW()
{
    return &W[RNG.irn(0, W.HowMany() - 1)];
}

vertex Diagram::RandomPickVer()
{
    return &Ver[RNG.irn(0, Ver.HowMany() - 1)];
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
        gLine g = G(index);
        vertex vin = NeighVer(g, IN);
        vertex vout = NeighVer(g, OUT);

        vin->nG[OUT] = g;
        vout->nG[IN] = g;

        g->Weight = GWeight->Weight(vin->R, vout->R, vin->Tau, vout->Tau, vin->Spin[OUT], vout->Spin[IN]);
        Weight *= g->Weight;
    }

    for (int index = 0; index < W.HowMany(); index++) {
        wLine w = W(index);
        vertex vin = NeighVer(w, IN);
        vertex vout = NeighVer(w, OUT);

        w->IsWorm = false;

        vin->nW = w;
        vout->nW = w;

        w->Weight = WWeight->Weight(vin->R, vout->R, vin->Tau, vout->Tau, vin->Spin, vout->Spin, w->IsWorm);
        Weight *= w->Weight;
    }

    for (int index = 0; index < Ver.HowMany(); index++) {
        //TODO: Do something here if you want to fix vertex
    }

    Phase = phase(Weight);
    return true;
}
