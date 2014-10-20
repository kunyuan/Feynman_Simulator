//
//  observable.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/13/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "weight.h"

using namespace std;
using namespace Array;
using namespace Weight;

Base::Base(const Lattice &lat, real beta)
{
    _Lat = lat;
    _Beta = beta;
}

int Base::SpinIndex(spin SpinIn, spin SpinOut)
{
    return SpinIn * SPIN + SpinOut;
}

int Base::SpinIndex(spin *TwoSpinIn, spin *TwoSpinOut)
{
    return TwoSpinIn[0] * SPIN3 + TwoSpinIn[1] * SPIN2 +
           TwoSpinOut[0] * SPIN + TwoSpinOut[1];
}

int Base::TauToBin(real tau)
{
    //TODO: mapping between tau and bin
    return int(tau * MAX_BIN / _Beta);
}

real Base::BinToTau(int Bin)
{
    //TODO: mapping between tau and bin
    return (real(Bin / MAX_BIN) + 0.5) * _Beta;
}

Sigma::Sigma(const Lattice &lat, real beta)
    : Base(lat, beta)
{
    _Shape[0] = SPIN2;
    _Shape[1] = lat.SublatVol * lat.SublatVol;
    _Shape[2] = lat.Vol;
    _Shape[3] = MAX_BIN;

    _Weight = new Array4<Complex>(_Shape[0], _Shape[1], _Shape[2], _Shape[3]);
    _WeightAccu = new Array4<Complex>(_Shape[0], _Shape[1], _Shape[2], _Shape[3]);
    _WeightSquareAccu = new Array4<Complex>(_Shape[0], _Shape[1], _Shape[2], _Shape[3]);

    _Name = "Sigma";
    _Norm = 1.0;
}

Sigma::~Sigma()
{
    delete _Weight;
    delete _WeightAccu;
    delete _WeightSquareAccu;
}

void Sigma::UpdateWeight()
{
    for (int i = 0; i < _Weight->Size(); i++) {
        (*_Weight)(i) = (*_WeightAccu)(i) / _Norm;
    }
}

Complex Sigma::Weight(const Distance &d, real dtau, spin SpinIn, spin SpinOut)
{
    return (*_Weight)[SpinIndex(SpinIn, SpinOut)][d.SublatIndex][d.CoordiIndex][TauToBin(dtau)];
}

Estimate<Complex> Sigma::WeightWithError(const Distance &d, real dtau, spin SpinIn, spin SpinOut)
{
    Complex sq2 = (*_WeightSquareAccu)[SpinIndex(SpinIn, SpinOut)][d.SublatIndex][d.CoordiIndex][TauToBin(dtau)] / _Norm;
    Complex mean = (*_WeightAccu)[SpinIndex(SpinIn, SpinOut)][d.SublatIndex][d.CoordiIndex][TauToBin(dtau)] / _Norm;
    return Estimate<Complex>(mean, sq2 - mean * mean);
}

void Sigma::Measure(const Complex &weight, const Distance &d, real dtau, spin SpinIn, spin SpinOut)
{
    int spin_index = SpinIndex(SpinIn, SpinOut);
    int tau_bin = TauToBin(dtau);
    (*_WeightAccu)[spin_index][d.SublatIndex][d.CoordiIndex][tau_bin] += weight;
    (*_WeightSquareAccu)[spin_index][d.SublatIndex][d.CoordiIndex][tau_bin] += weight * weight;
}

/************************   Polarization   *********************************/

Pi::Pi(const Lattice &lat, real beta)
    : Base(lat, beta)
{
    _Weight = new Array4<Complex>(SPIN4, lat.SublatVol * lat.SublatVol, lat.Vol, MAX_BIN);
    _WeightAccu = new Array4<Complex>(SPIN4, lat.SublatVol * lat.SublatVol,
                                      lat.Vol, MAX_BIN);
    _WeightSquareAccu = new Array4<Complex>(SPIN4, lat.SublatVol * lat.SublatVol,
                                            lat.Vol, MAX_BIN);
}

Pi::~Pi()
{
    delete _Weight;
    delete _WeightAccu;
    delete _WeightSquareAccu;
}

G::G(const Lattice &lat, real beta)
    : Base(lat, beta)
{
    _Weight = new Array4<Complex>(SPIN2, lat.SublatVol * lat.SublatVol, lat.Vol, MAX_BIN);
}

G::~G()
{
    delete _Weight;
}

Complex G::Weight(const Distance &dR, real dtau, spin SpinIn, spin SpinOut)
{
    return (*_Weight)[SpinIndex(SpinIn, SpinOut)][dR.SublatIndex][dR.CoordiIndex][TauToBin(dtau)];
}

W::W(const Lattice &lat, real beta)
    : Base(lat, beta)
{
    _Weight = new Array4<Complex>(SPIN4, lat.SublatVol * lat.SublatVol, lat.Vol, MAX_BIN);
}

Weight::W::~W()
{
    delete _Weight;
}
