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

Complex Sigma::Weight(const Site &rin, const Site& rout, real tin, real tout, spin SpinIn, spin SpinOut)
{
    return _Weight[SpinIndex(SpinIn, SpinOut)][_Lat.Dist(rin, rout).SublatIndex][_Lat.Dist(rin, rout).CoordiIndex][TauToBin(tout-tin)];
}

Complex Sigma::WeightOfDelta(spin SpinIn, spin SpinOut)
{
    //TODO: add Delta expression of Sigma here, you may need to add new API to measure it
    return 0.0;
}

void Sigma::Measure( const Site &rin, const Site& rout, real tin, real tout,spin SpinIn, spin SpinOut, int order, const Complex &weight)
{
    if (DEBUGMODE && order <= 0)
        LOG_ERROR("Too small order=" << order);
    int spin_index = SpinIndex(SpinIn, SpinOut);
    int tau_bin = TauToBin(tout-tin);
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

Complex Polar::Weight(const Site &rin, const Site& rout, real tin, real tout, spin *SpinIn, spin *SpinOut)
{
    return _Weight[SpinIndex(SpinIn, SpinOut)][_Lat.Dist(rin, rout).SublatIndex][_Lat.Dist(rin, rout).CoordiIndex][TauToBin(tout-tin)];
}

void Polar::Measure(const Site &rin, const Site& rout, real tin, real tout, spin *SpinIn, spin *SpinOut, int order, const Complex &weight)
{
    if (DEBUGMODE && order <= 0)
        LOG_ERROR("Too small order=" << order);
    int spin_index = SpinIndex(SpinIn, SpinOut);
    int tau_bin = TauToBin(tout-tin);
    _WeightAccu[order - 1][spin_index][_Lat.Dist(rin, rout).SublatIndex][_Lat.Dist(rin, rout).CoordiIndex][tau_bin] += weight;
    _Norm += _dBeta;
    _Average[order - 1].Measure(weight);
}

/***********************  G  **************************************/

G::G(const Lattice &lat, real beta, int order)
    : Weight::WeightNoMeasure(lat, beta, order, SPIN4, "Polar")
{
}

Complex G::Weight(const Site &rin, const Site& rout, real tin, real tout, spin SpinIn, spin SpinOut, bool IsMeasure)
{
    if(IsMeasure)
        //TODO: define the measuring weight for G
        return Complex(1.0, 0.0);
    else
        return _Weight[SpinIndex(SpinIn, SpinOut)][_Lat.Dist(rin, rout).SublatIndex][_Lat.Dist(rin, rout).CoordiIndex][TauToBin(tout-tin)];
}

Complex G::Weight(int dir, const Site &r1, const Site& r2, real t1, real t2, spin Spin1, spin Spin2, bool IsMeasure)
{
    if(IsMeasure)
        //TODO: define the measuring weight for G
        return Complex(1.0, 0.0);
    else
        if(dir==IN)
            return _Weight[SpinIndex(Spin1, Spin2)][_Lat.Dist(r1, r2).SublatIndex][_Lat.Dist(r1, r2).CoordiIndex][TauToBin(t2-t1)];
        else
            return _Weight[SpinIndex(Spin2, Spin1)][_Lat.Dist(r2, r1).SublatIndex][_Lat.Dist(r2, r1).CoordiIndex][TauToBin(t1-t2)];
}

Complex G::BareWeight(const Site &rin, const Site& rout, real tin, real tout, spin SpinIn, spin SpinOut)
{
    //TODO: add bare G expression here
    return 0.0;
}

/***********************  W  **************************************/

W::W(const Lattice &lat, real beta, int order)
    : WeightNoMeasure(lat, beta, order, SPIN4, "Polar")
{
}

Complex W::Weight(const Site &rin, const Site& rout, real tin, real tout, spin *SpinIn, spin *SpinOut, bool IsWorm, bool IsMeasure)
{
    if(IsMeasure)
        //TODO: define the measuring weight of W
        return Complex(1.0, 0.0);
    else
        if (IsWorm)
            //define your fake function here
            return _Weight[SpinIndex(UP, UP)][_Lat.Dist(rin, rout).SublatIndex][_Lat.Dist(rin, rout).CoordiIndex][TauToBin(tout-tin)];
        else
            return _Weight[SpinIndex(SpinIn, SpinOut)][_Lat.Dist(rin, rout).SublatIndex][_Lat.Dist(rin, rout).CoordiIndex][TauToBin(tout-tin)];
}

Complex W::Weight(int dir, const Site &r1, const Site& r2, real t1, real t2, spin *Spin1, spin *Spin2, bool IsWorm, bool IsMeasure)
{
    int spinindex, subindex, coordindex, tau;
    if(dir==IN)
    {
        spinindex = SpinIndex(Spin1, Spin2);
        subindex = _Lat.Dist(r1, r2).SublatIndex;
        coordindex = _Lat.Dist(r1, r2).CoordiIndex;
        tau = TauToBin(t2-t1);
    }else{
        spinindex = SpinIndex(Spin2, Spin1);
        subindex = _Lat.Dist(r2, r1).SublatIndex;
        coordindex = _Lat.Dist(r2, r1).CoordiIndex;
        tau = TauToBin(t1-t2);
    }
    
    
    if(IsMeasure)
        //TODO: define the measuring weight of W
        return Complex(1.0, 0.0);
    else
        if (IsWorm)
            //define your fake function here
            return _Weight[SpinIndex(UP, UP)][subindex][coordindex][tau];
        else
            return _Weight[spinindex][subindex][coordindex][tau];
}

Complex W::WeightOfDelta(const Site &rin, const Site& rout, spin *SpinIn, spin *SpinOut, bool IsWorm)
{
    //TODO: add Delta expression of W here, you may need to add new API to measure it
    return Complex(1.0,0.0);
}

Complex W::WeightOfDelta(int dir, const Site &r1, const Site& r2, spin *Spin1, spin *Spin2, bool IsWorm)
{
    //TODO: add Delta expression of W here, you may need to add new API to measure it
    return Complex(1.0,0.0);
}

Complex BareWeight(const Site &rin, const Site& rout, real tin, real tout, spin *SpinIn, spin *SpinOut)
{
    //TODO: add bare W expression here
    return 0.0;
}