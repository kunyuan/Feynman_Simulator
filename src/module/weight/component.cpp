//
//  component.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/22/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "component.h"
#include "utility/dictionary.h"
#include <tuple>

using namespace weight;
using namespace std;

real Norm::NormFactor = 1.0;

GClass::GClass(const Lattice& lat, real beta, uint MaxTauBin, TauSymmetry Symmetry)
    : _Map(IndexMapSPIN2(beta, MaxTauBin, lat, Symmetry))
{
    _SmoothTWeight.Allocate(_Map.GetShape(), SMOOTH);
    _SmoothTWeight.Assign(Complex(0.0, 0.0));
    _MeasureWeight.Allocate(_Map.GetShape(), SMOOTH);
    //initialize _MeasureWeight to an unit function
    _MeasureWeight.Assign(Complex(1.0, 0.0));
}

void GClass::BuildTest()
{
    _SmoothTWeight.Assign(0.0);
    for (uint sub = 0; sub < _Map.Lat.SublatVol; sub++) {
        Site Local(sub, { 0, 0 });
        for (uint tau = 0; tau < _Map.MaxTauBin; tau++) {
            Complex weight = exp(Complex(0.0, _Map.IndexToTau(tau)));
            uint Index = _Map.GetIndex(UP, UP, Local, Local, 0, tau);
            _SmoothTWeight[Index] = weight;
        }
    }
}

void GClass::Reset(real Beta)
{
    _Map = IndexMapSPIN2(Beta, _Map.MaxTauBin, _Map.Lat, _Map.Symmetry);
}

bool GClass::FromDict(const Dictionary& dict)
{
    return _SmoothTWeight.FromDict(dict);
}

Dictionary GClass::ToDict()
{
    return _SmoothTWeight.ToDict();
}

WClass::WClass(const Lattice& lat, real Beta, uint MaxTauBin)
    : _Map(IndexMapSPIN4(Beta, MaxTauBin, lat, TauSymmetric))
{
    _SmoothTWeight.Allocate(_Map.GetShape(), SMOOTH);
    _SmoothTWeight.Assign(Complex(0.0, 0.0));
    _DeltaTWeight.Allocate(_Map.GetShape(), DELTA);
    _DeltaTWeight.Assign(Complex(0.0, 0.0));
    _MeasureWeight.Allocate(_Map.GetShape(), SMOOTH);
    //initialize _MeasureWeight to an unit function
    _MeasureWeight.Assign(Complex(1.0, 0.0));
}

void WClass::BuildTest()
{
    _DeltaTWeight.Assign(0.0);
    _SmoothTWeight.Assign(0.0);
    spin UPUP[2] = { UP, UP };
    for (uint sub = 0; sub < _Map.Lat.SublatVol; sub++) {
        Site Local(sub, { 0, 0 });
        for (uint tau = 0; tau < _Map.MaxTauBin; tau++) {
            Complex weight = exp(Complex(0.0, -_Map.IndexToTau(tau)));
            uint Index = _Map.GetIndex(UPUP, UPUP, Local, Local, 0, tau);
            _SmoothTWeight[Index] = weight;
        }
    }
}

void WClass::Reset(real Beta)
{
    _Map = IndexMapSPIN4(Beta, _Map.MaxTauBin, _Map.Lat, _Map.Symmetry);
}

bool WClass::FromDict(const Dictionary& dict)
{
    return _SmoothTWeight.FromDict(dict) && _DeltaTWeight.FromDict(dict);
}

Dictionary WClass::ToDict()
{
    auto dict = _SmoothTWeight.ToDict();
    dict.Update(_DeltaTWeight.ToDict());
    return dict;
}

SigmaClass::SigmaClass(const Lattice& lat, real Beta, uint MaxTauBin,
             int MaxOrder, TauSymmetry Symmetry, real Norm)
    : _Map(IndexMapSPIN2(Beta, MaxTauBin, lat, Symmetry))
{
    Estimator.Allocate(_Map, MaxOrder, Norm);
}

void SigmaClass::BuildNew()
{
    Estimator.ClearStatistics();
}

void SigmaClass::BuildTest()
{
    Estimator.ClearStatistics();
}

void SigmaClass::Reset(real Beta)
{
    _Map = IndexMapSPIN2(Beta, _Map.MaxTauBin, _Map.Lat, _Map.Symmetry);
    Estimator.Anneal(Beta);
}

bool SigmaClass::FromDict(const Dictionary& dict)
{
    return Estimator.FromDict(dict.Get<Dictionary>("Histogram").Get<Dictionary>("SmoothT"));
}

Dictionary SigmaClass::ToDict()
{
    Dictionary dict;
    dict["Histogram"] = Dictionary("SmoothT", Estimator.ToDict());
    return dict;
}

PolarClass::PolarClass(const Lattice& lat, real Beta, uint MaxTauBin, int MaxOrder, real Norm)
    : _Map(IndexMapSPIN4(Beta, MaxTauBin, lat, TauSymmetric))
{
    Estimator.Allocate(_Map, MaxOrder, Norm);
}

void PolarClass::BuildNew()
{
    Estimator.ClearStatistics();
}

void PolarClass::BuildTest()
{
    Estimator.ClearStatistics();
}

void PolarClass::Reset(real Beta)
{
    _Map = IndexMapSPIN4(Beta, _Map.MaxTauBin, _Map.Lat, _Map.Symmetry);
    Estimator.Anneal(Beta);
}

bool PolarClass::FromDict(const Dictionary& dict)
{
    return Estimator.FromDict(dict.Get<Dictionary>("Histogram").Get<Dictionary>("SmoothT"));
}

Dictionary PolarClass::ToDict()
{
    Dictionary dict;
    dict["Histogram"] = Dictionary("SmoothT", Estimator.ToDict());
    return dict;
}
