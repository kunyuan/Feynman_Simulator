//
//  component.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/3/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "component.h"
#include "assert.h"
using namespace std;
using namespace diag;

const string ToString(const spin &s)
{
    if (s == UP)
        return "UP";
    else
        return "DOWN";
}

spin GLine::Spin(int dir)
{
    return nVer[dir]->_spin[INVERSE(dir)];
}

spin GLine::Spin()
{
    return nVer[0]->_spin[INVERSE(0)];
}

void GLine::FlipSpin()
{
    nVer[0]->_spin[1] = FLIP(nVer[0]->_spin[1]);
    nVer[1]->_spin[0] = FLIP(nVer[1]->_spin[0]);
}

int GLine::Sublattice(int dir)
{
    return nVer[dir]->R.Sublattice;
}

vertex GLine::NeighVer(int dir)
{
    return nVer[dir];
}

void GLine::SetGLine(Momentum k, const Complex& weight, bool ismeasure)
{
    K = k;
    Weight = weight;
    IsMeasure = ismeasure;
}

string GLine::PrettyString()
{
    stringstream os;
    os << "{V " << nVer[IN]->Name << "}->-" << ToString(Spin()) << "---";
    os << "[G " << Name << " ,K:" << K << ",Weight:" << Weight << "]";
    os << "---" << ToString(Spin()) << "->-{V " << nVer[OUT]->Name << "}";
    return os.str();
}

/*************** WLine  *************************/
spin WLine::Spin(int dir1, int dir2)
{
    return nVer[dir1]->_spin[dir2];
}

int WLine::Sublattice(int dir)
{
    return nVer[dir]->R.Sublattice;
}

vertex WLine::NeighVer(int dir)
{
    return nVer[dir];
}

void WLine::SetWLine(Momentum k, const Complex& weight, bool isworm, bool ismeasure, bool isdelta)
{
    K = k;
    Weight = weight;
    IsWorm = isworm;
    IsMeasure = ismeasure;
    IsDelta = isdelta;
}

string WLine::PrettyString()
{
    stringstream os;
    os << "{V " << nVer[IN]->Name << "| " << ToString(Spin(IN, IN)) << "," << ToString(Spin(IN, OUT)) << "}~~~";
    os << "<W " << Name << ", K:" << K << ",Weight:" << Weight << ">";
    os << "~~~"
       << "{" << ToString(Spin(OUT, IN)) << "," << ToString(Spin(OUT, OUT)) << "|V " << nVer[OUT]->Name << "}";
    return os.str();
}

/****************  Vertex **************************/

spin Vertex::Spin(int dir)
{
    return _spin[dir];
}

spin *Vertex::Spin()
{
    return _spin;
}

void Vertex::SetSpin(spin *_spin_)
{
    _spin[0] = _spin_[0];
    _spin[1] = _spin_[1];
}

int Vertex::Sublattice()
{
    return R.Sublattice;
}

gLine Vertex::NeighG(int dir)
{
    return nG[dir];
}

wLine Vertex::NeighW()
{
    return nW;
}

void Vertex::SetVertex(const Site& site, const real& tau, spin* s, int dir)
{
    R = site;
    Tau = tau;
    _spin[IN] = s[IN];
    _spin[OUT] = s[OUT];
    Dir = dir;
}

string Vertex::PrettyString()
{
    stringstream os;
    os << "-[G " << nG[IN]->Name << "]-" << ToString(Spin(IN)) << "->--";
    os << "{V " << Name << ",r:" << ToString(R.Coordinate) << ",tau:" << Tau << "}";
    os << "-->-" << ToString(Spin(OUT)) << "-[G " << nG[OUT]->Name << "]-";
    os << "  /~~~<W " << nW->Name << ">";
    return os.str();
}
