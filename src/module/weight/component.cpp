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
    //use _Shape[SP] to _Shape[TAU] to construct array3
    _MeasureWeight.Allocate(GetShape());
    //initialize _MeasureWeight to an unit function
    _MeasureWeight = Complex(1.0, 0.0);
}

void G::BuildTest(weight::model Model)
{
    if (Model == weight::Trivial) {
        _SmoothTWeight = 0.0;
        int spin_up = _Map.SpinIndex(UP, UP);
        for (int sub = 0; sub < GetShape()[SUB]; sub++) {
            if (!_Lat.IsOnSameSubLat(sub))
                continue;
            int coor = _Lat.Vec2Index({ 0, 0 });
            for (int tau = 0; tau < _Shape[TAU]; tau++) {
                Complex weight = exp(Complex(0.0, _Map.IndexToTau(tau)));
                _SmoothTWeight[spin_up][sub][coor][tau] = weight;
            }
        }
    }
    else if (Model == weight::DiagCount) {
        _SmoothTWeight = 0.0;
        int spin_up = _Map.SpinIndex(UP, UP);

        for (int sub = 0; sub < _Shape[SUB]; sub++) {
            if (!_Lat.IsOnSameSubLat(sub))
                continue;
            int coor = _Lat.Vec2Index({ 0, 0 });
            for (int tau = 0; tau < _Shape[TAU]; tau++) {
                Complex weight = Complex(1.0, 0.0);
                _SmoothTWeight[spin_up][sub][coor][tau] = weight;
            }
        }
    }
}

void G::Reset(real Beta)
{
    Basic::Reset(Beta);
    _Map = IndexMapSPIN2(Beta, _MaxTauBin, _Lat);
}

W::W(const Lattice& lat, real Beta, uint MaxTauBin)
    : weight::Basic(lat, Beta, MaxTauBin, SPIN4, TauSymmetric, "W")
    , _Map(IndexMapSPIN4(Beta, MaxTauBin, lat))
{
    //use _Shape[SP] to _Shape[VOL] to construct array3
    _MeasureWeight.Allocate(GetShape());
    //initialize _MeasureWeight to an unit function
    _MeasureWeight = Complex(1.0, 0.0);
}

void W::BuildTest(weight::model Model)
{
    if (Model == Trivial) {
        int Lx = _Lat.Size[0], Ly = _Lat.Size[1];
        ASSERT_ALLWAYS(Lx > 1 && Ly > 1, "System size should be bigger than 1!");
        int spin_up = _Map.SpinIndex(UP, //InOfW/InOfVertex
                                     UP, //InOfW/OutOfVertex
                                     UP, //OutOfW/InOfVertex
                                     UP); //OutOfW/OutOfVertex

        for (int sub = 0; sub < _Shape[SUB]; sub++) {
            if (!_Lat.IsOnSameSubLat(sub))
                continue;
            int coor = _Lat.Vec2Index({ 0, 0 });

            for (int tau = 0; tau < _Shape[TAU]; tau++) {
                Complex weight = exp(Complex(0.0, -_Map.IndexToTau(tau)));
                _SmoothTWeight[spin_up][sub][coor][tau] = weight;
            }
        }

        for (auto i : _Map.GetSpinIndexVector(UpUp2UpUp))
            _SmoothTWeight[i] = _SmoothTWeight[spin_up];

        for (auto i : _Map.GetSpinIndexVector(UpDown2UpDown)) {
            _SmoothTWeight[i] = _SmoothTWeight[spin_up];
            _SmoothTWeight[i] *= -1.0;
        }
        for (auto i : _Map.GetSpinIndexVector(UpDown2DownUp)) {
            _SmoothTWeight[i] = _SmoothTWeight[spin_up];
            _SmoothTWeight[i] *= 2.0;
        }
    }
    else if (Model == weight::DiagCount) {
        _DeltaTWeight = 0.0;
        _SmoothTWeight = 0.0;
        int Lx = _Lat.Size[0], Ly = _Lat.Size[1];
        ASSERT_ALLWAYS(Lx > 1 && Ly > 1, "System size should be bigger than 1!");
        int spin_up = _Map.SpinIndex(UP, //InOfW/InOfVertex
                                     UP, //InOfW/OutOfVertex
                                     UP, //OutOfW/InOfVertex
                                     UP); //OutOfW/OutOfVertex

        for (int sub = 0; sub < _Shape[SUB]; sub++) {
            if (!_Lat.IsOnSameSubLat(sub))
                continue;
            int coor = _Lat.Vec2Index({ 0, 0 });
            for (int tau = 0; tau < _Shape[TAU]; tau++) {
                Complex weight = Complex(1.0, 0.0);
                _SmoothTWeight[spin_up][sub][coor][tau] = weight;
            }
        }
    }
}

void W::Reset(real Beta)
{
    Basic::Reset(Beta);
    _Map = IndexMapSPIN4(Beta, _MaxTauBin, _Lat);
}

Sigma::Sigma(const Lattice& lat, real Beta, uint MaxTauBin,
             int MaxOrder, TauSymmetry Symmetry)
    : weight::Basic(lat, Beta, MaxTauBin, SPIN2, Symmetry, "Sigma")
    , _Map(IndexMapSPIN2(Beta, MaxTauBin, lat))
    , Estimator(lat, Beta, MaxOrder, "Sigma", Norm::Weight(), GetShape())
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
    Estimator.ReWeight(Beta);
}

bool Sigma::FromDict(const Dictionary& dict)
{
    return Estimator.FromDict(dict.Get<Dictionary>("Histogram").Get<Dictionary>("SmoothT"))
           && Basic::FromDict(dict);
}

Dictionary Sigma::ToDict()
{
    Dictionary dict = Basic::ToDict();
    dict["Histogram"] = Dictionary("SmoothT", Estimator.ToDict());
    return dict;
}

Polar::Polar(const Lattice& lat, real Beta, uint MaxTauBin, int MaxOrder)
    : weight::Basic(lat, Beta, MaxTauBin, SPIN4, TauSymmetric, "Polar")
    , _Map(IndexMapSPIN4(Beta, MaxTauBin, lat))
    , Estimator(lat, Beta, MaxOrder, "Polar", Norm::Weight(), GetShape())
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
    Estimator.ReWeight(Beta);
}

bool Polar::FromDict(const Dictionary& dict)
{
    return Estimator.FromDict(dict.Get<Dictionary>("Histogram").Get<Dictionary>("SmoothT"))
           && Basic::FromDict(dict);
}

Dictionary Polar::ToDict()
{
    Dictionary dict = Basic::ToDict();
    dict["Histogram"] = Dictionary("SmoothT", Estimator.ToDict());
    return dict;
}
