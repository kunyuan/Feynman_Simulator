//
//  component.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/11/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "diagram.h"
#include "../utility/abort.h"

#define SEP ' '
using namespace std;

#define READ(is, thing)                              \
    {                                                \
        is >> (thing);                               \
        if (is.fail())                               \
            ABORT("Fail to read " << #thing << "!"); \
    }

/*******************   Read/write component to dat file  ********************************/
bool Diagram::LoadConfig(istream &is, WormClass &worm)
{
    //format: i/Ira/Masha/dSpin/K
    name ira, masha;
    READ(is, ira);
    worm.Ira = Ver(ira);
    READ(is, masha);
    worm.Masha = Ver(masha);
    READ(is, worm.dSpin);
    READ(is, worm.K);
    return true;
}

void Diagram::SaveConfig(ostream &os, WormClass &worm)
{
    os << worm.Ira->Name << SEP << worm.Masha->Name << SEP << worm.dSpin << SEP << worm.K << endl;
}

bool Diagram::LoadConfig(istream &is, gLine g)
{
    //format: g/start/end/K
    name in, out;
    READ(is, in);
    g->nVer[IN] = Ver(in);
    READ(is, out);
    g->nVer[OUT] = Ver(out);
    READ(is, g->K);
    return true;
}

void Diagram::SaveConfig(ostream &os, gLine g)
{
    os << g->nVer[IN]->Name << SEP << g->nVer[OUT]->Name << SEP << g->K << endl;
}

bool Diagram::LoadConfig(istream &is, wLine w)
{
    //format: w/start/end
    name in, out;
    READ(is, in);
    w->nVer[IN] = Ver(in);
    READ(is, out);
    w->nVer[OUT] = Ver(out);
    READ(is, w->K);
    return true;
}

void Diagram::SaveConfig(ostream &os, wLine w)
{
    os << w->nVer[IN]->Name << SEP << w->nVer[OUT]->Name << SEP << w->K << endl;
}

bool Diagram::LoadConfig(istream &is, vertex v)
{
    int spinin, spinout;
    //format: name sublattice r tau
    READ(is, v->Name);
    READ(is, v->R.Sublattice);
    READ(is, v->R.Coordinate);
    READ(is, v->Tau);
    READ(is, spinin);
    READ(is, spinout);
    if (is.good()) {
        v->Spin[IN] = spin(spinin);
        v->Spin[OUT] = spin(spinout);
    }
    return true;
}

void Diagram::SaveConfig(ostream &os, vertex v)
{
    os << v->Name << SEP << v->R.Sublattice << SEP << v->R.Coordinate << SEP << v->Tau << SEP << int(v->Spin[IN]) << SEP << int(v->Spin[OUT]) << endl;
}
