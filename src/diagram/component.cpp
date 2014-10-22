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
ostream &operator<<(ostream &os, GLine &r)
{
    os << r.Vertex[IN] << SEP << r.Vertex[OUT] << SEP<< r.K << endl;
    return os;
}

istream &operator>>(istream &is, GLine &r)
{
    //format: start/end/in_spin/out_spin
    is >> r.Vertex[IN] >> r.Vertex[OUT] >>r.K;
    return is;
}
ostream &operator<<(ostream &os, WLine &r)
{
    os << r.Vertex[IN] << SEP << r.Vertex[OUT] <<SEP << r.K<< endl;
    return os;
}

istream &operator>>(istream &is, WLine &r)
{
    //format: start/end
    is >> r.Vertex[IN] >> r.Vertex[OUT] >>r.K;
    return is;
}

ostream &operator<<(ostream &os, Vertex &r)
{
    os << r.Name << SEP << r.R.Sublattice << SEP << r.R.Coordinate << SEP << r.Tau << SEP << int(r.Spin[IN]) << SEP << int(r.Spin[OUT]) << endl;
    return os;
}

istream &operator>>(istream &is, Vertex &r)
{
    int spinin, spinout;
    //format: name sublattice r tau
    is >> r.Name >> r.R.Sublattice >> r.R.Coordinate >> r.Tau >> spinin >> spinout;
    r.Spin[IN]=spin(spinin);
    r.Spin[OUT]=spin(spinout);
    return is;
}
