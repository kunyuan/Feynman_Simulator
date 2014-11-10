//
//  diagram_global.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/9/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

//TODO: G, W IsDelta, IsMeasure

#include "diagram.h"
#include "module/weight/weight.h"
#include "utility/rng.h"

using namespace std;
using namespace diag;

bool Diagram::IsWorm(vertex v)
{
    if (Worm.Ira == v || Worm.Masha == v)
        return true;
    else
        return false;
}

Diagram::Diagram()
    : Order(0), Phase(Complex(1.0, 0.0)), Weight(Complex(1.0, 0.0)), G("GLine"), W("WLine"), Ver("nVer")
{
    Lat = nullptr;
    GWeight = nullptr;
    WWeight = nullptr;
}

#include "diagram_initialize.config"
void Diagram::BuildNew(Lattice &lat, RandomFactory &rng, weight::G *g, weight::W *w)
{
    Reset(lat, rng, g, w);
    stringstream ss(InitialDiagram);
    _Load(ss);
    FixDiagram();
}

void Diagram::Reset(Lattice &lat, RandomFactory &rng, weight::G *g, weight::W *w)
{
    Lat = &lat;
    RNG = &rng;
    GWeight = g;
    WWeight = w;
    FixDiagram();
    //TODO: maybe you have to do more to reset
}

#include "diagram_template.config"
void Diagram::SetTest(Lattice &lat, RandomFactory &rng, weight::G *g, weight::W *w)
{
    Reset(lat, rng, g, w);
    stringstream ss(TestDiagramString);
    _Load(ss);
    FixDiagram();
}

gLine Diagram::RandomPickG()
{
    return &G[RNG->irn(0, G.HowMany() - 1)];
}

wLine Diagram::RandomPickW()
{
    return &W[RNG->irn(0, W.HowMany() - 1)];
}

vertex Diagram::RandomPickVer()
{
    return &Ver[RNG->irn(0, Ver.HowMany() - 1)];
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
    Worm.Weight = 1.0;

    Weight = Complex(1.0, 0.0);
    for (int index = 0; index < G.HowMany(); index++) {
        gLine g = G(index);
        vertex vin = g->NeighVer(IN);
        vertex vout = g->NeighVer(OUT);

        vin->nG[OUT] = g;
        vout->nG[IN] = g;

        g->Weight = GWeight->Weight(vin->R, vout->R, vin->Tau, vout->Tau, vin->Spin(OUT), vout->Spin(IN), g->IsMeasure);
        Weight *= g->Weight;
    }

    for (int index = 0; index < W.HowMany(); index++) {
        wLine w = W(index);
        vertex vin = w->NeighVer(IN);
        vertex vout = w->NeighVer(OUT);

        w->IsWorm = false;

        vin->nW = w;
        vin->Dir = IN;

        vout->nW = w;
        vout->Dir = OUT;

        w->Weight = WWeight->Weight(vin->R, vout->R, vin->Tau, vout->Tau, vin->Spin(), vout->Spin(), w->IsWorm, w->IsMeasure, w->IsDelta);
        Weight *= w->Weight;
    }

    for (int index = 0; index < Ver.HowMany(); index++) {
        //TODO: Do something here if you want to fix vertex
    }

    Phase = SignFermiLoop * (Order%2==0? 1:-1) *phase(Weight);
    return true;
}
