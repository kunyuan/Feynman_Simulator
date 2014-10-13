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
ostream &operator<<(ostream &os, spin &s)
{
    if (s == UP)
        return os << "UP";
    else
        return os << "DOWN";
}

string GLine::PrettyString()
{
    stringstream os;
    os << "\n";
    os << Vertex[IN] << " (" << Spin[IN] << ") >>===";
    os << "Name:" << Name << ",Weight:" << Weight;
    os << "===>>" << Vertex[OUT] << " (" << Spin[OUT] << ");";
    os << endl;
    return os.str();
}

ostream &operator<<(ostream &os, GLine &r)
{
    os << r.Vertex[IN] << SEP << r.Vertex[OUT] << SEP;
    os << int(r.Spin[IN]) << SEP << int(r.Spin[OUT]) << endl;
    return os;
}

istream &operator>>(istream &is, GLine &r)
{
    //format: start/end/in_spin/out_spin
    int spin_in, spin_out;
    is >> r.Vertex[IN] >> r.Vertex[OUT] >> spin_in >> spin_out;
    r.Spin[IN] = (spin)spin_in;
    r.Spin[OUT] = (spin)spin_out;
    return is;
}

string WLine::PrettyString()
{
    stringstream os;
    os << "\n";
    os << Vertex[IN] << " (" << Spin[IN][IN] << "," << Spin[IN][OUT] << ") ~~~";
    os << "Name:" << Name << ",Weight:" << Weight;
    os << "~~~" << Vertex[OUT] << " (" << Spin[OUT][IN] << "," << Spin[OUT][IN] << ");";
    os << endl;
    return os.str();
}

ostream &operator<<(ostream &os, WLine &r)
{
    os << r.Vertex[IN] << SEP << r.Vertex[OUT] << endl;
    return os;
}

istream &operator>>(istream &is, WLine &r)
{
    //format: start/end
    is >> r.Vertex[IN] >> r.Vertex[OUT];
    return is;
}

string Vertex::PrettyString()
{
    stringstream os;
    os << "\n";
    os << G[IN] << ">===";
    os << "Name:" << Name << ",r:" << r << ",tau:" << tau;
    os << "===>" << G[OUT] << "\n";
    os << "       ~~~~" << W;
    os << endl;
    return os.str();
}

ostream &operator<<(ostream &os, Vertex &r)
{
    os << r.Name << SEP << r.Sublattice << SEP << r.r << SEP << r.tau << endl;
    return os;
}

istream &operator>>(istream &is, Vertex &r)
{
    //format: name/sublattice/r/tau
    int cor;
    is >> r.Name >> r.Sublattice >> cor >> r.tau;
    return is;
}
