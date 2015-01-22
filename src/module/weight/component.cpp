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

G::G(const Lattice& lat, real beta, uint MaxTauBin, TauSymmetry Symmetry)
    : _Map(IndexMapSPIN2(beta, MaxTauBin, lat, Symmetry))
{
    _SmoothTWeight.Allocate(_Map.GetShape(), SMOOTH);
    _SmoothTWeight.Assign(Complex(0.0, 0.0));
    _MeasureWeight.Allocate(_Map.GetShape(), SMOOTH);
    //initialize _MeasureWeight to an unit function
    _MeasureWeight.Assign(Complex(1.0, 0.0));
}

void G::BuildTest()
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

void G::Reset(real Beta)
{
    _Map = IndexMapSPIN2(Beta, _Map.MaxTauBin, _Map.Lat, _Map.Symmetry);
}

bool G::FromDict(const Dictionary& dict)
{
    return _SmoothTWeight.FromDict(dict);
}

Dictionary G::ToDict()
{
    return _SmoothTWeight.ToDict();
}

W::W(const Lattice& lat, real Beta, uint MaxTauBin)
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

void W::BuildTest()
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

void W::Reset(real Beta)
{
    _Map = IndexMapSPIN4(Beta, _Map.MaxTauBin, _Map.Lat, _Map.Symmetry);
}

bool W::FromDict(const Dictionary& dict)
{
    return _SmoothTWeight.FromDict(dict) && _DeltaTWeight.FromDict(dict);
}

Dictionary W::ToDict()
{
    auto dict = _SmoothTWeight.ToDict();
    dict.Update(_DeltaTWeight.ToDict());
    return dict;
}

Sigma::Sigma(const Lattice& lat, real Beta, uint MaxTauBin,
             int MaxOrder, TauSymmetry Symmetry, real Norm)
    : _Map(IndexMapSPIN2(Beta, MaxTauBin, lat, Symmetry))
{
    Estimator.Allocate(_Map, MaxOrder, Norm);
}

void Sigma::BuildNew()
{
    Estimator.ClearStatistics();
}

void Sigma::BuildTest()
{
    Estimator.ClearStatistics();
}

void Sigma::Reset(real Beta)
{
    _Map = IndexMapSPIN2(Beta, _Map.MaxTauBin, _Map.Lat, _Map.Symmetry);
    Estimator.Anneal(Beta);
}

bool Sigma::FromDict(const Dictionary& dict)
{
    return Estimator.FromDict(dict.Get<Dictionary>("Histogram").Get<Dictionary>("SmoothT"));
}

Dictionary Sigma::ToDict()
{
    Dictionary dict;
    dict["Histogram"] = Dictionary("SmoothT", Estimator.ToDict());
    return dict;
}

Polar::Polar(const Lattice& lat, real Beta, uint MaxTauBin, int MaxOrder, real Norm)
    : _Map(IndexMapSPIN4(Beta, MaxTauBin, lat, TauSymmetric))
{
    Estimator.Allocate(_Map, MaxOrder, Norm);
}

void Polar::BuildNew()
{
    Estimator.ClearStatistics();
}

void Polar::BuildTest()
{
    Estimator.ClearStatistics();
}

void Polar::Reset(real Beta)
{
    _Map = IndexMapSPIN4(Beta, _Map.MaxTauBin, _Map.Lat, _Map.Symmetry);
    Estimator.Anneal(Beta);
}

bool Polar::FromDict(const Dictionary& dict)
{
    return Estimator.FromDict(dict.Get<Dictionary>("Histogram").Get<Dictionary>("SmoothT"));
}

Dictionary Polar::ToDict()
{
    Dictionary dict;
    dict["Histogram"] = Dictionary("SmoothT", Estimator.ToDict());
    return dict;
}
