//
//  weight_IO.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/20/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "weight.h"

using namespace std;
using namespace Array;
using namespace Weight;

Sigma::Sigma(const Lattice &lat, real beta, int order)
    : WeightNeedMeasure(lat, beta, order, SPIN2, "Sigma")
{
}

Complex Sigma::Weight(const Distance &d, real dtau, spin SpinIn, spin SpinOut)
{
    return _Weight[SpinIndex(SpinIn, SpinOut)][d.SublatIndex][d.CoordiIndex][TauToBin(dtau)];
}

Complex Sigma::WeightOfDelta(spin SpinIn, spin SpinOut)
{
    //TODO: add Delta expression of Sigma here, you may need to add new API to measure it
    return 0.0;
}

void Sigma::Measure(const Distance &d, real dtau, spin SpinIn, spin SpinOut, int order, const Complex &weight)
{
    if (DEBUGMODE && order <= 0)
        LOG_ERROR("Too small order=" << order);
    int spin_index = SpinIndex(SpinIn, SpinOut);
    int tau_bin = TauToBin(dtau);
    _WeightAccu[order - 1][spin_index][d.SublatIndex][d.CoordiIndex][tau_bin] += weight;
    _Norm += _dBeta;
    _Average[order - 1].Measure(weight);
}

/************************   Polarization   *********************************/
//
Polar::Polar(const Lattice &lat, real beta, int order)
    : WeightNeedMeasure(lat, beta, order, SPIN4, "Polar")
{
}

Complex Polar::Weight(const Distance &d, real dtau, spin *SpinIn, spin *SpinOut)
{
    return _Weight[SpinIndex(SpinIn, SpinOut)][d.SublatIndex][d.CoordiIndex][TauToBin(dtau)];
}

void Polar::Measure(const Distance &d, real dtau, spin *SpinIn, spin *SpinOut, int order, const Complex &weight)
{
    if (DEBUGMODE && order <= 0)
        LOG_ERROR("Too small order=" << order);
    int spin_index = SpinIndex(SpinIn, SpinOut);
    int tau_bin = TauToBin(dtau);
    _WeightAccu[order - 1][spin_index][d.SublatIndex][d.CoordiIndex][tau_bin] += weight;
    _Norm += _dBeta;
    _Average[order - 1].Measure(weight);
}

/***********************  G  **************************************/

G::G(const Lattice &lat, real beta, int order)
    : Weight::WeightNoMeasure(lat, beta, order, SPIN4, "Polar")
{
}

Complex G::Weight(const Distance &d, real dtau, spin SpinIn, spin SpinOut)
{
    return _Weight[SpinIndex(SpinIn, SpinOut)][d.SublatIndex][d.CoordiIndex][TauToBin(dtau)];
}

Complex G::BareWeight(const Distance &d, real dtau, spin SpinIn, spin SpinOut)
{
    //TODO: add bare G expression here
    return 0.0;
}

/***********************  W  **************************************/

W::W(const Lattice &lat, real beta, int order)
    : WeightNoMeasure(lat, beta, order, SPIN4, "Polar")
{
}

Complex W::Weight(const Distance &d, real dtau, spin *SpinIn, spin *SpinOut, bool IsWorm)
{
    if (IsWorm)
        //define your fake function here
        return _Weight[SpinIndex(UP, UP)][d.SublatIndex][d.CoordiIndex][TauToBin(dtau)];
    else
        return _Weight[SpinIndex(SpinIn, SpinOut)][d.SublatIndex][d.CoordiIndex][TauToBin(dtau)];
}

Complex W::WeightOfDelta(const Distance &d, spin *SpinIn, spin *SpinOut, bool IsWorm)
{
    //TODO: add Delta expression of W here, you may need to add new API to measure it
    return 0.0;
}

Complex BareWeight(const Distance &dR, real dtau, spin *SpinIn, spin *SpinOut)
{
    //TODO: add bare W expression here
    return 0.0;
}