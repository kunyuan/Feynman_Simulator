//
//  component_MC.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/24/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "component.h"

using namespace weight0;
using namespace std;

#define GetSmoothT(array, index) array[index[SP]][index[SUB]][index[VOL]][index[TAU]]
#define GetDeltaT(array, index) array[index[SP]][index[SUB]][index[TAU]]
const spin SPINUPUP[2] = {UP, UP};

Complex G::Weight(const Site &rin, const Site &rout, real tin, real tout, spin SpinIn, spin SpinOut, bool IsMeasure)
{
    static uint index[4];
    _Map.Map(index, SpinIn, SpinOut, rin, rout, tin, tout);
    if (IsMeasure)
        return GetSmoothT(_MeasureWeight, index);
    else
        return _TauSymmetryFactor * GetSmoothT(_SmoothTWeight, index);
}

Complex G::Weight(int dir, const Site &r1, const Site &r2, real t1, real t2, spin Spin1, spin Spin2, bool IsMeasure)
{
    static uint index[4];
    if (dir == IN)
        _Map.Map(index, Spin1, Spin2, r1, r2, t1, t2);
    else
        _Map.Map(index, Spin2, Spin1, r2, r1, t2, t1);

    if (IsMeasure)
        return GetSmoothT(_MeasureWeight, index);
    else
        return _TauSymmetryFactor * GetSmoothT(_SmoothTWeight, index);
}

Complex W::Weight(const Site &rin, const Site &rout, real tin, real tout, spin *SpinIn, spin *SpinOut, bool IsWorm, bool IsMeasure, bool IsDelta)
{
    static uint index[4];
    if (IsWorm) {
        //it is safe to reassign pointer here, the original spins pointed by SpinIn and SpinOut pointers will not change
        SpinIn = (spin *)SPINUPUP;
        SpinOut = (spin *)SPINUPUP;
    }
    if (IsDelta) {
        _Map.MapDeltaT(index, SpinIn, SpinOut, rin, rout);
        return GetDeltaT(_DeltaTWeight, index);
    }
    _Map.Map(index, SpinIn, SpinOut, rin, rout, tin, tout);
    if (IsMeasure)
        return GetSmoothT(_MeasureWeight, index);
    else
        return GetSmoothT(_SmoothTWeight, index);
}

Complex W::Weight(int dir, const Site &r1, const Site &r2, real t1, real t2, spin *Spin1, spin *Spin2, bool IsWorm, bool IsMeasure, bool IsDelta)
{
    static uint index[4];
    if (IsWorm) {
        Spin1 = (spin *)SPINUPUP;
        Spin2 = (spin *)SPINUPUP;
    }
    if (dir == IN)
        if (IsDelta)
            _Map.MapDeltaT(index, Spin1, Spin2, r1, r2);
        else
            _Map.Map(index, Spin1, Spin2, r1, r2, t1, t2);
    else if (IsDelta)
        _Map.MapDeltaT(index, Spin2, Spin1, r2, r1);
    else
        _Map.Map(index, Spin2, Spin1, r2, r1, t2, t1);

    if (IsMeasure)
        return GetSmoothT(_MeasureWeight, index);
    else if (IsDelta)
        return GetDeltaT(_DeltaTWeight, index);
    else
        return GetSmoothT(_SmoothTWeight, index);
}

Complex Sigma::Weight(const Site &rin, const Site &rout, real tin, real tout, spin SpinIn, spin SpinOut)
{
    static uint index[4];
    _Map.Map(index, SpinIn, SpinOut, rin, rout, tin, tout);
    return _TauSymmetryFactor * GetSmoothT(_SmoothTWeight, index);
}

void Sigma::Measure(const Site &rin, const Site &rout, real tin, real tout, spin SpinIn, spin SpinOut, int order, const Complex &weight)
{
    static uint index[4];
    _Map.Map(index, SpinIn, SpinOut, rin, rout, tin, tout);
    Estimator.Measure(index, order, weight * GetTauSymmetryFactor(tin, tout));
}

Complex Polar::Weight(const Site &rin, const Site &rout, real tin, real tout, spin *SpinIn, spin *SpinOut)
{
    static uint index[4];
    _Map.Map(index, SpinIn, SpinOut, rin, rout, tin, tout);
    return GetSmoothT(_SmoothTWeight, index);
}

void Polar::Measure(const Site &rin, const Site &rout, real tin, real tout, spin *SpinIn, spin *SpinOut, int order, const Complex &weight)
{
    static uint index[4];
    _Map.Map(index, SpinIn, SpinOut, rin, rout, tin, tout);
    Estimator.Measure(index, order, weight);
}