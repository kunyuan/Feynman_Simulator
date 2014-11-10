//
//  weight_IO.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/20/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "weight_inherit.h"

using namespace std;
using namespace Array;
using namespace weight;

Sigma::Sigma(const Lattice &lat, real beta, int order)
    : WeightNeedMeasure(lat, beta, order, SPIN2, "Sigma", Norm::Weight())
{
}

Complex Sigma::Weight(const Site &rin, const Site &rout, real tin, real tout, spin SpinIn, spin SpinOut)
{
    return SmoothWeight[SpinIndex(SpinIn, SpinOut)][_Lat.Dist(rin, rout).SublatIndex][_Lat.Dist(rin, rout).CoordiIndex][TauToBin(tout - tin)];
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
    int spin = SpinIndex(SpinIn, SpinOut);
    auto dist = _Lat.Dist(rin, rout);
    int tau = TauToBin(tin, tout);
    _WeightAccu[order - 1][spin][dist.SublatIndex][dist.CoordiIndex][tau] += weight;
    if (spin == 0 && dist.SublatIndex == 0 && dist.CoordiIndex == 0 && tau == 0)
        _WeightErrorEstimator[order - 1].Measure(weight);
}

/************************   Polarization   *********************************/
//
Polar::Polar(const Lattice &lat, real beta, int order)
    : WeightNeedMeasure(lat, beta, order, SPIN4, "Polar", Norm::Weight())
{
}

Complex Polar::Weight(const Site &rin, const Site &rout, real tin, real tout, spin *SpinIn, spin *SpinOut)
{
    return SmoothWeight[SpinIndex(SpinIn, SpinOut)][_Lat.Dist(rin, rout).SublatIndex][_Lat.Dist(rin, rout).CoordiIndex][TauToBin(tout - tin)];
}

void Polar::Measure(const Site &rin, const Site &rout, real tin, real tout, spin *SpinIn, spin *SpinOut, int order, const Complex &weight)
{
    if (DEBUGMODE && order <= 0)
        LOG_ERROR("Too small order=" << order);
    int spin = SpinIndex(SpinIn, SpinOut);
    auto dist = _Lat.Dist(rin, rout);
    int tau = TauToBin(tin, tout);
    _WeightAccu[order - 1][spin][dist.SublatIndex][dist.CoordiIndex][tau] += weight;
    if (spin == 0 && dist.SublatIndex == 0 && dist.CoordiIndex == 0 && tau == 0)
        _WeightErrorEstimator[order - 1].Measure(weight);
}

/***********************  G  **************************************/

G::G(const Lattice &lat, real beta, int order)
    : WeightNoMeasure(lat, beta, order, SPIN4, "G")
{
}

Complex G::Weight(const Site &rin, const Site &rout, real tin, real tout, spin SpinIn, spin SpinOut, bool IsMeasure)
{
    auto dist = _Lat.Dist(rin, rout);
    if (IsMeasure)
        //TODO: define the measuring weight for G
        return Complex(1.0, 0.0);
    else
        return SmoothWeight[SpinIndex(SpinIn, SpinOut)]
                           [dist.SublatIndex]
                           [dist.CoordiIndex]
                           [TauToBin(tin, tout)];
}

Complex G::Weight(int dir, const Site &r1, const Site &r2, real t1, real t2, spin Spin1, spin Spin2, bool IsMeasure)
{
    if (IsMeasure)
        //TODO: define the measuring weight for G
        return Complex(1.0, 0.0);

    else if (dir == IN) {
        auto dist = _Lat.Dist(r1, r2);
        return SmoothWeight[SpinIndex(Spin1, Spin2)]
                           [dist.SublatIndex]
                           [dist.CoordiIndex]
                           [TauToBin(t1, t2)];
    }
    else {
        auto dist = _Lat.Dist(r2, r1);
        return SmoothWeight[SpinIndex(Spin2, Spin1)]
                           [dist.SublatIndex]
                           [dist.CoordiIndex]
                           [TauToBin(t2, t1)];
    }
}

void G::InitialWithBare()
{
    //TODO: add bare G initialization
}

/***********************  W  **************************************/

W::W(const Lattice &lat, real beta, int order)
    : WeightNoMeasure(lat, beta, order, SPIN4, "W")
{
}

Complex W::Weight(const Site &rin, const Site &rout, real tin, real tout, spin *SpinIn, spin *SpinOut, bool IsWorm, bool IsMeasure, bool IsDelta)
{
    if (IsMeasure)
        //TODO: define the measuring weight of W
        return Complex(1.0, 0.0);

    auto distance = _Lat.Dist(rin, rout);
    if (IsDelta)
        //TODO: define the delta function here! IsWorm==true and IsWorm==false
        return DeltaTWeight[SpinIndex(SpinIn, SpinOut)]
                           [distance.SublatIndex]
                           [distance.CoordiIndex];

    else if (IsWorm)
        //define your fake function here
        return SmoothWeight[SpinIndex(UP, UP)]
                           [distance.SublatIndex]
                           [distance.CoordiIndex]
                           [TauToBin(tout - tin)];
    else
        return SmoothWeight[SpinIndex(SpinIn, SpinOut)]
                           [distance.SublatIndex]
                           [distance.CoordiIndex]
                           [TauToBin(tin, tout)];
}

Complex W::Weight(int dir, const Site &r1, const Site &r2, real t1, real t2, spin *Spin1, spin *Spin2, bool IsWorm, bool IsMeasure, bool IsDelta)
{
    int spinindex, subindex, coordindex, tau;
    if (dir == IN) {
        auto dist = _Lat.Dist(r1, r2);
        spinindex = SpinIndex(Spin1, Spin2);
        subindex = dist.SublatIndex;
        coordindex = dist.CoordiIndex;
        tau = TauToBin(t1, t2);
    }
    else {
        auto dist = _Lat.Dist(r2, r1);
        spinindex = SpinIndex(Spin2, Spin1);
        subindex = dist.SublatIndex;
        coordindex = dist.CoordiIndex;
        tau = TauToBin(t2, t1);
    }
    if (IsMeasure)
        //TODO: define the measuring weight of W
        return Complex(1.0, 0.0);
    if (IsDelta)
        //TODO: define the delta function here! IsWorm==true and IsWorm==false
        return Complex(1.0, 0.0);
    else if (IsWorm)
        //define your fake function here
        return SmoothWeight[SpinIndex(UP, UP)][subindex][coordindex][tau];
    else
        return SmoothWeight[spinindex][subindex][coordindex][tau];
}

void W::InitialWithBare()
{
    DeltaTWeight = 1.0;
    SmoothWeight = 0.0;
    //TODO: add bare W initialization
}