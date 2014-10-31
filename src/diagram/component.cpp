//
//  component.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/11/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "component.h"
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
istream &Worm::LoadConfig(istream &is)
{
    //format: i/Ira/Masha/dSpin/K
    READ(is, Ira->Name);
    READ(is, Masha);
    READ(is, dSpin);
    READ(is, K);
    return is;
}

ostream &Worm::SaveConfig(ostream &os)
{
    os << Ira << SEP << Masha << SEP << dSpin << SEP << K << endl;
    return os;
}

istream &GLine::LoadConfig(istream &is)
{
    //format: g/start/end/in_spin/out_spin
    size_t vertex_in,vertex_out;
    READ(is, vertex_in);
    nVer[IN]=(Vertex *)vertex_in;
    READ(is, vertex_out);
    nVer[OUT]=(Vertex *)vertex_out;
    READ(is, K);
    return is;
}

ostream &GLine::SaveConfig(ostream &os)
{
    os << nVer[IN] << SEP << nVer[OUT] << SEP << K << endl;
    return os;
}

istream &WLine::LoadConfig(istream &is)
{
    //format: start/end
    READ(is, Vertex[IN]);
    READ(is, Vertex[OUT]);
    READ(is, K);
    return is;
}

ostream &WLine::SaveConfig(ostream &os)
{
    os << Vertex[IN] << SEP << Vertex[OUT] << SEP << K << endl;
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
