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
    : weight::Basic(lat, beta, MaxTauBin, SPIN2, Symmetry, "G")
    , _Map(IndexMapSPIN2(beta, MaxTauBin, lat))
{
    _SmoothTWeight.Allocate(GetShape());
    _SmoothTWeight.Assign(Complex(0.0, 0.0));
    //use _Shape[SP] to _Shape[TAU] to construct array3
    _MeasureWeight.Allocate(GetShape());
    //initialize _MeasureWeight to an unit function
    _MeasureWeight.Assign(Complex(1.0, 0.0));
}

void G::BuildTest()
{
    _SmoothTWeight.Assign(0.0);
    uint spin_up = _Map.SpinIndex(UP, UP);
    for (uint sub = 0; sub < GetShape()[SUB]; sub++) {
        if (!_Lat.IsOnSameSubLat(sub))
            continue;
        uint coor = _Lat.Vec2Index({ 0, 0 });
        for (uint tau = 0; tau < _Shape[TAU]; tau++) {
            Complex weight = exp(Complex(0.0, _Map.IndexToTau(tau)));
            _SmoothTWeight({ spin_up, sub, coor, tau }) = weight;
        }
    }
}

void G::Reset(real Beta)
{
    Basic::Reset(Beta);
    _Map = IndexMapSPIN2(Beta, _MaxTauBin, _Lat);
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
    : weight::Basic(lat, Beta, MaxTauBin, SPIN4, TauSymmetric, "W")
    , _Map(IndexMapSPIN4(Beta, MaxTauBin, lat))
{
    _SmoothTWeight.Allocate(GetShape());
    _SmoothTWeight.Assign(Complex(0.0, 0.0));
    _DeltaTWeight.Allocate(GetShape());
    _DeltaTWeight.Assign(Complex(0.0, 0.0));
    //use _Shape[SP] to _Shape[VOL] to construct array3
    _MeasureWeight.Allocate(GetShape());
    //initialize _MeasureWeight to an unit function
    _MeasureWeight.Assign(Complex(1.0, 0.0));
}

void W::BuildTest()
{
    int Lx = _Lat.Size[0], Ly = _Lat.Size[1];
    ASSERT_ALLWAYS(Lx > 1 && Ly > 1, "System size should be bigger than 1!");
    uint spin_up = _Map.SpinIndex(UP, //InOfW/InOfVertex
                                  UP, //InOfW/OutOfVertex
                                  UP, //OutOfW/InOfVertex
                                  UP); //OutOfW/OutOfVertex

    for (uint sub = 0; sub < _Shape[SUB]; sub++) {
        if (!_Lat.IsOnSameSubLat(sub))
            continue;
        uint coor = _Lat.Vec2Index({ 0, 0 });

        for (uint tau = 0; tau < _Shape[TAU]; tau++) {
            Complex weight = exp(Complex(0.0, -_Map.IndexToTau(tau)));
            _SmoothTWeight({ spin_up, sub, coor, tau }) = weight;
        }
    }
}

void W::Reset(real Beta)
{
    Basic::Reset(Beta);
    _Map = IndexMapSPIN4(Beta, _MaxTauBin, _Lat);
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
    : weight::Basic(lat, Beta, MaxTauBin, SPIN2, Symmetry, "Sigma")
    , _Map(IndexMapSPIN2(Beta, MaxTauBin, lat))
    , Estimator(lat, Beta, MaxOrder, "Sigma", Norm, GetShape())
{
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
    Basic::Reset(Beta);
    _Map = IndexMapSPIN2(Beta, _MaxTauBin, _Lat);
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
    : weight::Basic(lat, Beta, MaxTauBin, SPIN4, TauSymmetric, "Polar")
    , _Map(IndexMapSPIN4(Beta, MaxTauBin, lat))
    , Estimator(lat, Beta, MaxOrder, "Polar", Norm, GetShape())
{
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
    Basic::Reset(Beta);
    _Map = IndexMapSPIN4(Beta, _MaxTauBin, _Lat);
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
