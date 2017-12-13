//
//  diagram_global.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/9/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

//TODO: G, W IsDelta, IsMeasure

#include "diagram.h"
#include "module/weight/component.h"
#include "utility/rng.h"

using namespace std;
using namespace diag;
using namespace weight;

bool Diagram::IsWorm(vertex v)
{
    if (Worm.Exist && (Worm.Ira == v || Worm.Masha == v))
        return true;
    else
        return false;
}

Diagram::Diagram()
    : Order(0)
    , Phase(Complex(1.0, 0.0))
    , Weight(Complex(1.0, 0.0))
    , G("GLine")
    , W("WLine")
    , Ver("nVer")
    , HasGammaGW(0)
{
    Lat = nullptr;
}

bool Diagram::GHashCheck(Momentum k)
{
    return (GHash[k.index()]);
}

void Diagram::AddGHash(Momentum k)
{
    if (DEBUGMODE && GHash[k.index()])
        ABORT("add occupied G Hash!");
    GHash[k.index()] = true;
}

void Diagram::RemoveGHash(Momentum k)
{
    if (DEBUGMODE && !(GHash[k.index()]))
        ABORT("remove empty G Hash!");
    GHash[k.index()] = false;
}

void Diagram::ReplaceGHash(Momentum kold, Momentum k)
{
    RemoveGHash(kold);
    AddGHash(k);
}

bool Diagram::WHashCheck(Momentum k)
{
    return (k == 0 || WHash[k.abs()]);
}

void Diagram::AddWHash(Momentum k)
{
    if (DEBUGMODE && WHash[k.abs()])
        ABORT("add occupied W Hash!");
    WHash[k.abs()] = true;
}

void Diagram::RemoveWHash(Momentum k)
{
    if (DEBUGMODE && !(WHash[k.abs()]))
        ABORT("remove empty W Hash!");
    WHash[k.abs()] = false;
}

void Diagram::ReplaceWHash(Momentum kold, Momentum k)
{
    RemoveWHash(kold);
    AddWHash(k);
}

void Diagram::Reset(Lattice& lat, weight::GClass& g, weight::WClass& w)
{
    Lat = &lat;
    GWeight = &g;
    WWeight = &w;
    LOG_WARNING("Calling fix diagram!");
    FixDiagram();
}

void Diagram::ClearDiagram()
{
    for (int i = 0; i < 2 * MAX_K + 1; i++) {
        GHash[i] = false;
    }
    for (int i = 0; i < MAX_K + 1; i++) {
        WHash[i] = false;
    }

    while (G.HowMany() > 0)
        G.Remove(G.HowMany() - 1);
    while (W.HowMany() > 0)
        W.Remove(W.HowMany() - 1);
    while (Ver.HowMany() > 0)
        Ver.Remove(Ver.HowMany() - 1);
}
bool Diagram::FixDiagram()
{

    if(HasGammaGW!=0)
        ABORT("Fix diagram when gammaG and GammaW exists!");

    ASSERT_ALLWAYS(Lat != nullptr, "Lattice is not defined yet!");
    ASSERT_ALLWAYS(GWeight != nullptr && WWeight != nullptr, "G and W have been initialized yet!");

    Order = W.HowMany();
    Worm.Weight = 1.0;
    
    if(Worm.Exist){
        Worm.Weight = weight::Worm::Weight(Worm.Ira->R, Worm.Masha->R, Worm.Ira->Tau, Worm.Masha->Tau);
    }

    Weight = Complex(1.0, 0.0);
    for (int index = 0; index < G.HowMany(); index++) {
        gLine g = G(index);
        vertex vin = g->NeighVer(IN);
        vertex vout = g->NeighVer(OUT);

        vin->nG[OUT] = g;
        vout->nG[IN] = g;
        g->IsGGGammaG = false;

        g->Weight = GWeight->Weight(vin->R, vout->R, vin->Tau, vout->Tau, vin->Spin(OUT), vout->Spin(IN), g->IsMeasure, false);
        Weight *= g->Weight;
    }

    for (int index = 0; index < W.HowMany(); index++) {
        wLine w = W(index);
        vertex vin = w->NeighVer(IN);
        vertex vout = w->NeighVer(OUT);

        w->IsWorm = (vin==Worm.Ira || vin==Worm.Masha || vout==Worm.Ira ||vout==Worm.Masha);

        vin->nW = w;
        vin->Dir = IN;

        vout->nW = w;
        vout->Dir = OUT;
        w->IsWWGammaW = false;

        w->Weight = WWeight->Weight(vin->R, vout->R, vin->Tau, vout->Tau, vin->Spin(), vout->Spin(), w->IsWorm, w->IsMeasure, w->IsDelta, false);
        Weight *= w->Weight;
    }

    for (int index = 0; index < Ver.HowMany(); index++) {
        //TODO: Do something here if you want to fix vertex
    }

    Weight *= SignFermiLoop * (Order % 2 == 0 ? 1 : -1);
    Phase = phase(Weight);
    return true;
}
