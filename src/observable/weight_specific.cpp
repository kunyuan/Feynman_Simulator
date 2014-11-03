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
using namespace weight;

Sigma::Sigma(const Lattice &lat, real beta, int order)
    : WeightNeedMeasure(lat, beta, order, SPIN2, "Sigma")
{
}

Complex Sigma::Weight(const Site &rin, const Site &rout, real tin, real tout, spin SpinIn, spin SpinOut)
{
    return _Weight[SpinIndex(SpinIn, SpinOut)][_Lat.Dist(rin, rout).SublatIndex][_Lat.Dist(rin, rout).CoordiIndex][TauToBin(tout - tin)];
}

Complex Sigma::WeightOfDelta(spin SpinIn, spin SpinOut)
{
    //TODO: add Delta expression of Sigma here, you may need to add new API to measure it
    return 0.0;
}

void Sigma::Measure(const Site &rin, const Site &rout, real tin, real tout, spin SpinIn, spin SpinOut, int order, const Complex &weight)
{
    if (DEBUGMODE && order <= 0)
        LOG_ERROR("Too small order=" << order);
    int spin_index = SpinIndex(SpinIn, SpinOut);
    int tau_bin = TauToBin(tout - tin);
    _WeightAccu[order - 1][spin_index][_Lat.Dist(rin, rout).SublatIndex][_Lat.Dist(rin, rout).CoordiIndex][tau_bin] += weight;
    _Norm += _dBeta;
    _Average[order - 1].Measure(weight);
}

/************************   Polarization   *********************************/
//
Polar::Polar(const Lattice &lat, real beta, int order)
    : WeightNeedMeasure(lat, beta, order, SPIN4, "Polar")
{
}

Complex Polar::Weight(const Site &rin, const Site &rout, real tin, real tout, spin *SpinIn, spin *SpinOut)
{
    return _Weight[SpinIndex(SpinIn, SpinOut)][_Lat.Dist(rin, rout).SublatIndex][_Lat.Dist(rin, rout).CoordiIndex][TauToBin(tout - tin)];
}

void Polar::Measure(const Site &rin, const Site &rout, real tin, real tout, spin *SpinIn, spin *SpinOut, int order, const Complex &weight)
{
    if (DEBUGMODE && order <= 0)
        LOG_ERROR("Too small order=" << order);
    int spin_index = SpinIndex(SpinIn, SpinOut);
    int tau_bin = TauToBin(tout - tin);
    _WeightAccu[order - 1][spin_index][_Lat.Dist(rin, rout).SublatIndex][_Lat.Dist(rin, rout).CoordiIndex][tau_bin] += weight;
    _Norm += _dBeta;
    _Average[order - 1].Measure(weight);
}

/***********************  G  **************************************/

G::G(const Lattice &lat, real beta, int order)
    : WeightNoMeasure(lat, beta, order, SPIN4, "Polar")
{
}

Complex G::Weight(const Site &rin, const Site &rout, real tin, real tout, spin SpinIn, spin SpinOut)
{
    return _Weight[SpinIndex(SpinIn, SpinOut)][_Lat.Dist(rin, rout).SublatIndex][_Lat.Dist(rin, rout).CoordiIndex][TauToBin(tout - tin)];
}

void G::InitialWithBare()
{
    //TODO: add bare G expression here
}

/***********************  W  **************************************/

W::W(const Lattice &lat, real beta, int order)
    : WeightNoMeasure(lat, beta, order, SPIN4, "Polar")
{
}

Complex W::Weight(const Site &rin, const Site &rout, real tin, real tout, spin *SpinIn, spin *SpinOut, bool IsWorm)
{
    if (IsWorm)
        //define your fake function here
        return _Weight[SpinIndex(UP, UP)][_Lat.Dist(rin, rout).SublatIndex][_Lat.Dist(rin, rout).CoordiIndex][TauToBin(tout - tin)];
    else
        return _Weight[SpinIndex(SpinIn, SpinOut)][_Lat.Dist(rin, rout).SublatIndex][_Lat.Dist(rin, rout).CoordiIndex][TauToBin(tout - tin)];
}

Complex W::WeightOfDelta(const Site &rin, const Site &rout, spin *SpinIn, spin *SpinOut, bool IsWorm)
{
    //TODO: add Delta expression of W here, you may need to add new API to measure it
    return 0.0;
}

void W::InitialWithBare()
{
    //TODO: add bare W expression here
}