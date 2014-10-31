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
    READ(is, worm.Ira);
    READ(is, worm.Masha);
    READ(is, worm.dSpin);
    READ(is, worm.K);
    return true;
}

void Diagram::SaveConfig(ostream &os, WormClass &worm)
{
    os << worm.Ira << SEP << worm.Masha << SEP << worm.dSpin << SEP << worm.K << endl;
}

bool Diagram::LoadConfig(istream &is, gLine g)
{
    //format: g/start/end/in_spin/out_spin
    size_t in, out;
    READ(is, in);
    g->nVer[IN] = Ver.ToPointer(in);
    READ(is, out);
    g->nVer[OUT] = Ver.ToPointer(out);
    READ(is, g->K);
    return true;
}

void Diagram::SaveConfig(ostream &os, gLine g)
{
    os << Ver.ToIndex(g->nVer[IN]) << SEP << Ver.ToIndex(g->nVer[OUT]) << SEP << g->K << endl;
}

istream &WLine::LoadConfig(istream &is)
{
    //format: start/end
    READ(is, nVer[IN]);
    READ(is, nVer[OUT]);
    READ(is, K);
    return is;
}

ostream &WLine::SaveConfig(ostream &os)
{
    os << nVer[IN] << SEP << nVer[OUT] << SEP << K << endl;
    return os;
}

istream &Vertex::LoadConfig(istream &is)
{
    int spinin, spinout;
    //format: name sublattice r tau
    READ(is, Name);
    READ(is, R.Sublattice);
    READ(is, R.Coordinate);
    READ(is, Tau);
    READ(is, spinin);
    READ(is, spinout);
    if (is.good()) {
        Spin[IN] = spin(spinin);
        Spin[OUT] = spin(spinout);
    }
    return is;
}

ostream &Vertex::SaveConfig(ostream &os)
{
    os << Name << SEP << R.Sublattice << SEP << R.Coordinate << SEP << Tau << SEP << int(Spin[IN]) << SEP << int(Spin[OUT]) << endl;
    return os;
}
