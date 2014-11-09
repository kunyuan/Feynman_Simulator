//
//  weightMonteCarlo.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/9/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "weightMonteCarlo.h"

using namespace std;
using namespace Array;
using namespace weight0;

mc::Sigma::Sigma(const Lattice &lat, real beta, int order)
    : WeightEstimator(lat, beta, SPIN2, "Sigma", order,
                      Norm::Weight() //set WeightEstimator._Norm
                      )
{
}

void mc::Sigma::Measure(const Site &rin, const Site &rout, real tin, real tout, spin SpinIn, spin SpinOut, int order, const Complex &weight)
{
    if (DEBUGMODE && order <= 0)
        LOG_ERROR("Too small order=" << order);
    auto dist = _Lat.Dist(rin, rout);
    _WeightAccu[order - 1]
               [SpinIndex(SpinIn, SpinOut)]
               [dist.SublatIndex]
               [dist.CoordiIndex]
               [TauToBin(tin, tout)] += weight;
    _NormAccu += 1.0;
    _Average[order - 1].Measure(weight);
}

/************************   Polarization   *********************************/
//
mc::Polar::Polar(const Lattice &lat, real beta, int order)
    : WeightEstimator(lat, beta, SPIN4, "Polar", order,
                      Norm::Weight() //set WeightEstimator._Norm
                      )
{
}

void mc::Polar::Measure(const Site &rin, const Site &rout, real tin, real tout, spin *SpinIn, spin *SpinOut, int order, const Complex &weight)
{
    if (DEBUGMODE && order <= 0)
        LOG_ERROR("Too small order=" << order);
    auto dist = _Lat.Dist(rin, rout);
    _WeightAccu[order - 1]
               [SpinIndex(SpinIn, SpinOut)]
               [dist.SublatIndex]
               [dist.CoordiIndex]
               [TauToBin(tin, tout)] += weight;
    _NormAccu += 1.0;
    _Average[order - 1].Measure(weight);
}

/***********************  G  **************************************/

mc::G::G(const Lattice &lat, real beta)
    : WeightArray(lat, beta, SPIN4, "G")
{
}

Complex mc::G::Weight(const Site &rin, const Site &rout, real tin, real tout, spin SpinIn, spin SpinOut, bool IsMeasure)
{
    auto dist = _Lat.Dist(rin, rout);
    if (IsMeasure)
        //TODO: define the measuring weight for G
        return Complex(1.0, 0.0);
    else
        return WeightArray::Weight[SpinIndex(SpinIn, SpinOut)]
                                  [dist.SublatIndex]
                                  [dist.CoordiIndex]
                                  [TauToBin(tin, tout)];
}

Complex mc::G::Weight(int dir, const Site &r1, const Site &r2, real t1, real t2, spin Spin1, spin Spin2, bool IsMeasure)
{
    if (IsMeasure)
        //TODO: define the measuring weight for G
        return Complex(1.0, 0.0);

    else if (dir == IN) {
        auto dist = _Lat.Dist(r1, r2);
        return WeightArray::Weight[SpinIndex(Spin1, Spin2)]
                                  [dist.SublatIndex]
                                  [dist.CoordiIndex]
                                  [TauToBin(t1, t2)];
    }
    else {
        auto dist = _Lat.Dist(r2, r1);
        return WeightArray::Weight[SpinIndex(Spin2, Spin1)]
                                  [dist.SublatIndex]
                                  [dist.CoordiIndex]
                                  [TauToBin(t2, t1)];
    }
}

void mc::G::InitialWithBare()
{
    //TODO: add bare G initialization
}

/***********************  W  **************************************/

mc::W::W(const Lattice &lat, real beta)
    : WeightArray(lat, beta, SPIN4, "W")
{
}

Complex mc::W::Weight(const Site &rin, const Site &rout, real tin, real tout, spin *SpinIn, spin *SpinOut, bool IsWorm, bool IsMeasure, bool IsDelta)
{
    if (IsMeasure)
        //TODO: define the measuring weight of W
        return Complex(1.0, 0.0);

    auto distance = _Lat.Dist(rin, rout);
    if (IsDelta)
        //TODO: define the delta function here! IsWorm==true and IsWorm==false
        return WeightArray::DeltaTWeight[SpinIndex(SpinIn, SpinOut)]
                                        [distance.SublatIndex]
                                        [distance.CoordiIndex];

    else if (IsWorm)
        //define your fake function here
        return WeightArray::Weight[SpinIndex(UP, UP)]
                                  [distance.SublatIndex]
                                  [distance.CoordiIndex]
                                  [TauToBin(tout - tin)];
    else
        return WeightArray::Weight[SpinIndex(SpinIn, SpinOut)]
                                  [distance.SublatIndex]
                                  [distance.CoordiIndex]
                                  [TauToBin(tin, tout)];
}

Complex mc::W::Weight(int dir, const Site &r1, const Site &r2, real t1, real t2, spin *Spin1, spin *Spin2, bool IsWorm, bool IsMeasure, bool IsDelta)
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
        return WeightArray::Weight[SpinIndex(UP, UP)][subindex][coordindex][tau];
    else
        return WeightArray::Weight[spinindex][subindex][coordindex][tau];
}

void mc::W::InitialWithBare()
{
    WeightArray::DeltaTWeight = 1.0;
    WeightArray::Weight = 0.0;
    //TODO: add bare W initialization
}