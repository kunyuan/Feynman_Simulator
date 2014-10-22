//
//  component.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/11/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "component.h"
#include <sstream>

#define SEP ' '
using namespace std;

/*******************   Read/write component to dat file  ********************************/
istream &Worm::LoadConfig(istream &is)
{
    //format: i/Ira/Masha/dSpin/K
    is >> Ira >> Masha >> dSpin >> K;
    return is;
}

ostream &Worm::SaveConfig(ostream &os)
{
    os << Ira <<SEP<< Masha <<SEP<< dSpin << SEP<<K <<endl;
    return os;
}

istream &GLine::LoadConfig(istream &is)
{
    //format: g/start/end/in_spin/out_spin
    is >> Vertex[IN] >> Vertex[OUT] >> K;
    return is;
}

ostream &GLine::SaveConfig(ostream &os)
{
    os << Vertex[IN] << SEP << Vertex[OUT] << SEP << K << endl;
    return os;
}

istream &WLine::LoadConfig(istream &is)
{
    //format: start/end
    is >> Vertex[IN] >> Vertex[OUT] >> K;
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
    is >> Name >> R.Sublattice >> R.Coordinate >> Tau >> spinin >> spinout;
    Spin[IN] = spin(spinin);
    Spin[OUT] = spin(spinout);
    return is;
}

ostream &Vertex::SaveConfig(ostream &os)
{
    os << Name << SEP << R.Sublattice << SEP << R.Coordinate << SEP << Tau << SEP << int(Spin[IN]) << SEP << int(Spin[OUT]) << endl;
    return os;
}
