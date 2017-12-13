//
//  component_MC.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/24/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "component.h"
#include "gamma3.h"

using namespace weight;
using namespace std;

const spin SPINUPUP[2] = { UP, UP };

Complex GClass::Weight(const Site& rin, const Site& rout, real tin, real tout, spin SpinIn, spin SpinOut,
                       bool IsMeasure, bool IsGGGammaG, ExtPoint* ext) const
{
    uint Index = _Map.GetIndex(SpinIn, SpinOut, rin, rout, tin, tout);
    if (IsMeasure)
        return _MeasureWeight(Index);
    else if (IsGGGammaG){
        return GGGammaGWeight->Weight(rin, rout, ext->R, tin, tout, ext->Tau, SpinIn, SpinOut, ext->Spin);
//        if(ext->R == rin)
//            return _Map.GetTauSymmetryFactor(tin, tout) * _SmoothTWeight(Index);
//        else
//            return Complex(0.0, 0.0);
    }
    else
        return _Map.GetTauSymmetryFactor(tin, tout) * _SmoothTWeight(Index);
}

Complex GClass::Weight(int dir, const Site& r1, const Site& r2, real t1, real t2, spin Spin1, spin Spin2,
                       bool IsMeasure, bool IsGGGammaG, ExtPoint* ext) const
{
    static uint Index;
    int symmetryfactor;
    if (dir == IN) {
        Index = _Map.GetIndex(Spin1, Spin2, r1, r2, t1, t2);
        symmetryfactor = _Map.GetTauSymmetryFactor(t1, t2);
    }
    else {
        Index = _Map.GetIndex(Spin2, Spin1, r2, r1, t2, t1);
        symmetryfactor = _Map.GetTauSymmetryFactor(t2, t1);
    }

    if (IsMeasure)
        return _MeasureWeight(Index);
    else if (IsGGGammaG)
//        if(ext->R == r1)
//            return symmetryfactor * _SmoothTWeight(Index);
//        else
//            return Complex(0.0, 0.0);
        if (dir == IN) {
            return GGGammaGWeight->Weight(r1, r2, ext->R, t1, t2, ext->Tau, Spin1, Spin2, ext->Spin);
        }else{
            return GGGammaGWeight->Weight(r2, r1, ext->R, t2, t1, ext->Tau, Spin2, Spin1, ext->Spin);
        }
    else
        return symmetryfactor * _SmoothTWeight(Index);
}

Complex WClass::Weight(const Site& rin, const Site& rout, real tin, real tout, spin* SpinIn, spin* SpinOut, bool IsWorm,
                       bool IsMeasure, bool IsDelta, bool IsWWGammaW, ExtPoint* ext) const
{
    static uint index;
    if (IsWorm) {
        //it is safe to reassign pointer here, the original spins pointed by SpinIn and SpinOut pointers will not change
        SpinIn = (spin*)SPINUPUP;
        SpinOut = (spin*)SPINUPUP;
    }
    if (IsDelta) {
        index = _Map.GetIndex(SpinIn, SpinOut, rin, rout);
        return _DeltaTWeight(index);
    }
    index = _Map.GetIndex(SpinIn, SpinOut, rin, rout, tin, tout);

    if (IsMeasure)
        return _MeasureWeight(index);
    else if(IsWWGammaW)
        return WWGammaWWeight->Weight(rin, rout, ext->R, tin, tout, ext->Tau, SpinIn, SpinOut, ext->Spin);
    else
        return _SmoothTWeight(index);
}

Complex WClass::Weight(int dir, const Site& r1, const Site& r2, real t1, real t2, spin* Spin1, spin* Spin2, bool IsWorm,
                       bool IsMeasure, bool IsDelta, bool IsWWGammaW, ExtPoint* ext) const
{
    static uint index;
    if (IsWorm) {
        Spin1 = (spin*)SPINUPUP;
        Spin2 = (spin*)SPINUPUP;
    }
    if (dir == IN)
        if (IsDelta)
            index = _Map.GetIndex(Spin1, Spin2, r1, r2);
        else
            index = _Map.GetIndex(Spin1, Spin2, r1, r2, t1, t2);

    else if (IsDelta)
        index = _Map.GetIndex(Spin2, Spin1, r2, r1);
    else
        index = _Map.GetIndex(Spin2, Spin1, r2, r1, t2, t1);

    if (IsMeasure)
        return _MeasureWeight(index);
    else if (IsDelta)
        return _DeltaTWeight(index);
    else if(IsWWGammaW)
        if (dir == IN)
            return WWGammaWWeight->Weight(r1, r2, ext->R, t1, t2, ext->Tau, Spin1, Spin2, ext->Spin);
        else
            return WWGammaWWeight->Weight(r2, r1, ext->R, t2, t1, ext->Tau, Spin2, Spin1, ext->Spin);
    else
        return _SmoothTWeight(index);
}

void SigmaClass::Measure(const Site& rin, const Site& rout, real tin, real tout, spin SpinIn, spin SpinOut, int order, const Complex& weight)
{
    static uint index;
    index = _Map.GetIndex(SpinIn, SpinOut, rin, rout, tin, tout);
    Estimator.Measure(index, order, weight * _Map.GetTauSymmetryFactor(tin, tout));
}

void PolarClass::Measure(const Site& rin, const Site& rout, real tin, real tout, spin* SpinIn, spin* SpinOut, int order, const Complex& weight)
{
    static uint index;
    index = _Map.GetIndex(SpinIn, SpinOut, rin, rout, tin, tout);
    Estimator.Measure(index, order, weight);
}
