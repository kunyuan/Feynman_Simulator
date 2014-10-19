//
//  observable.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/13/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "diagram_object.h"
using namespace Array;

inline int SpinIndex(spin SpinIn, spin SpinOut)
{
    return SpinIn*SPIN+SpinOut;
}

inline int SpinIndex(spin* TwoSpinIn, spin* TwoSpinOut)
{
    return TwoSpinIn[0]*SPIN3+TwoSpinIn[1]*SPIN2+
        TwoSpinOut[0]*SPIN+TwoSpinOut[1];
}

inline int TauToBin_Sym(real tau)
{
    return 0;
}

SigmaWeight::SigmaWeight(int Vol,int Sublattice2)
{
    _Weight=new Array4<Complex>(SPIN2,Sublattice2,Vol,MAX_BIN);
    _WeightAccu=new Array4<Complex>(SPIN2,Sublattice2,Vol,MAX_BIN);
    _WeightSquareAccu=new Array4<Complex>(SPIN2,Sublattice2,Vol,MAX_BIN);
}

SigmaWeight::~SigmaWeight()
{
    delete[] _Weight;
    delete[] _WeightAccu;
    delete[] _WeightSquareAccu;
}

void SigmaWeight::UpdateWeight()
{
    for(int i=0;i<_Weight->Size();i++)
    {
        (*_Weight)(i)=(*_WeightAccu)(i)/_Norm;
    }
}

Complex SigmaWeight::Weight(const Distance& dR, real dtau, spin SpinIn, spin SpinOut)
{
    return (*_Weight)[SpinIndex(SpinIn, SpinOut)][dR.Sublattice()][dR.Coordinate()][TauToBin_Sym(dtau)];
}

Estimate<Complex> SigmaWeight::WeightWithError(const Distance& dR, real dtau, spin SpinIn, spin SpinOut)
{
    Complex sq2=(*_WeightSquareAccu)[SpinIndex(SpinIn, SpinOut)][dR.Sublattice()][dR.Coordinate()][TauToBin_Sym(dtau)]/_Norm;
    Complex mean= (*_WeightAccu)[SpinIndex(SpinIn, SpinOut)][dR.Sublattice()][dR.Coordinate()][TauToBin_Sym(dtau)]/_Norm;
    return Estimate<Complex>(mean, sq2-mean*mean);
}

void SigmaWeight::Measure(const Complex &weight,const Distance& dR, real dtau, spin SpinIn, spin SpinOut)
{
    int spin_index=SpinIndex(SpinIn, SpinOut);
    int tau_bin=TauToBin_Sym(dtau);
    _WeightAccu[spin_index][dR.Sublattice()][dR.Coordinate()][tau_bin]+=weight;
    _WeightSquareAccu[spin_index][dR.Sublattice()][dR.Coordinate()][tau_bin]+=weight*weight;
}

