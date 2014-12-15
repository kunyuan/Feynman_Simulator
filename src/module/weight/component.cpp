//
//  component.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/22/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "component.h"
#include "utility/cnpy.h"
#include <tuple>

using namespace weight;
using namespace std;

G::G(const Lattice& lat, real beta, TauSymmetry Symmetry)
    : weight::Basic(lat, beta, SPIN2, Symmetry, "G")
    , _Map(IndexMapSPIN2(beta, lat))
{
    //use _Shape[SP] to _Shape[TAU] to construct array3
    _MeasureWeight.Allocate(GetShape());
    //initialize _MeasureWeight to an unit function
    _MeasureWeight = Complex(1.0, 0.0);
}

void G::BuildNew(model Model,
                 const std::vector<real>& Hopping,
                 const std::vector<Complex>& ChemicalPotential,
                 real ExternalField)
{
    _Hopping = Hopping;
    _ExternalField = ExternalField;
    _ChemicalPotential = ChemicalPotential;
    _Model = Model;
    _SmoothTWeight = 0.0;
    int spin_down = _Map.SpinIndex(DOWN, DOWN);
    int spin_up = _Map.SpinIndex(UP, UP);

    for (int sub = 0; sub < _Shape[SUB]; sub++) {
        if (!_Lat.IsOnSameSubLat(sub))
            continue;
        int coor = _Lat.Vec2Index({ 0, 0 });
        //        MeasureWeight[spin_up][sub][coor] = Complex(1.0, 0.0);
        for (int tau = 0; tau < _Shape[TAU]; tau++) {
            Complex weight = Complex(1.0, 0.0);
            _SmoothTWeight[spin_down][sub][coor][tau] = weight;
            _SmoothTWeight[spin_up][sub][coor][tau] = weight;
        }
    }
}

void G::BuildTest()
{
    _Model = model::TEST;
    _SmoothTWeight = 0.0;
    int spin_down = _Map.SpinIndex(DOWN, DOWN);
    int spin_up = _Map.SpinIndex(UP, UP);
    for (int sub = 0; sub < GetShape()[SUB]; sub++) {
        if (!_Lat.IsOnSameSubLat(sub))
            continue;
        int coor = _Lat.Vec2Index({ 0, 0 });
        for (int tau = 0; tau < _Shape[TAU]; tau++) {
            Complex weight = exp(Complex(0.0, _Map.IndexToTau(tau)));
            _SmoothTWeight[spin_down][sub][coor][tau] = weight;
            _SmoothTWeight[spin_up][sub][coor][tau] = weight;
        }
    }
}

void G::Reset(real Beta)
{
    Basic::Reset(Beta);
    _Map = IndexMapSPIN2(Beta, _Lat);
}

W::W(const Lattice& lat, real Beta)
    : weight::Basic(lat, Beta, SPIN4, TauSymmetric, "W")
    , _Map(IndexMapSPIN4(Beta, lat))
{
    //use _Shape[SP] to _Shape[VOL] to construct array3
    _MeasureWeight.Allocate(GetShape());
    //initialize _MeasureWeight to an unit function
    _MeasureWeight = Complex(1.0, 0.0);
}

void W::BuildNew(model Model, const vector<real>& Interaction_, real ExternalField)
{
    _Interaction = Interaction_;
    _ExternalField = ExternalField;
    _Model = Model;
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

void W::BuildTest()
{
    _Model = model::TEST;
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

void W::Reset(real Beta)
{
    Basic::Reset(Beta);
    _Map = IndexMapSPIN4(Beta, _Lat);
}

void W::WriteBareToASCII()
{
    _Lat.PlotLattice();
    Vec<int> offset;
    for (int i = 0; i < D; i++)
        offset[i] = _Lat.Size[i] / 2 - 1;
    auto Shape = GetShape();
    ofstream os("interaction.py", ios::out);
    int spin_index = _Map.SpinIndex(UP, UP, UP, UP);
    os << "line=[" << endl;
    for (int sub = 0; sub < Shape[SUB]; sub++)
        for (int coord = 0; coord < Shape[VOL]; coord++) {
            real bare_weight = mod(_DeltaTWeight[spin_index][sub][coord]);
            if (!Zero(bare_weight)) {
                Site start, end;
                tie(start, end) = _Lat.GetSite(Distance(sub, coord));
                os << "[" << _Lat.GetRealVec(start, offset).PrettyString() << ","
                   << _Lat.GetRealVec(end, offset).PrettyString() << ","
                   << start.Sublattice << "]," << endl;
            }
        }
    os << "]" << endl;
}

Sigma::Sigma(const Lattice& lat, real Beta, int MaxOrder, TauSymmetry Symmetry)
    : weight::Basic(lat, Beta, SPIN2, Symmetry, "Sigma")
    , _Map(IndexMapSPIN2(Beta, lat))
    , Estimator(Beta, MaxOrder, "Sigma", Norm::Weight(), GetShape())
{
}

void Sigma::BuildNew(model Model)
{
    _Model = Model;
    Estimator.ClearStatistics();
}

void Sigma::BuildTest()
{
    _Model = model::TEST;
    Estimator.ClearStatistics();
}

void Sigma::Reset(real Beta)
{
    Basic::Reset(Beta);
    _Map = IndexMapSPIN2(Beta, _Lat);
    Estimator.ReWeight(Beta);
}

bool Sigma::Load(const std::string& FileName)
{
    return Estimator.Load(FileName) && Basic::Load(FileName);
}
void Sigma::Save(const std::string& FileName, const std::string Mode)
{
    Estimator.Save(FileName, Mode);
    Basic::Save(FileName, "a");
}

Polar::Polar(const Lattice& lat, real Beta, int MaxOrder)
    : weight::Basic(lat, Beta, SPIN4, TauSymmetric, "Polar")
    , _Map(IndexMapSPIN4(Beta, lat))
    , Estimator(Beta, MaxOrder, "Polar", Norm::Weight(), GetShape())
{
}

void Polar::BuildNew(model Model)
{
    _Model = Model;
    Estimator.ClearStatistics();
}

void Polar::BuildTest()
{
    _Model = model::TEST;
    Estimator.ClearStatistics();
}

void Polar::Reset(real Beta)
{
    Basic::Reset(Beta);
    _Map = IndexMapSPIN4(Beta, _Lat);
    Estimator.ReWeight(Beta);
}

bool Polar::Load(const std::string& FileName)
{
    return Estimator.Load(FileName) && Basic::Load(FileName);
}
void Polar::Save(const std::string& FileName, const std::string Mode)
{
    Estimator.Save(FileName, Mode);
    Basic::Save(FileName, "a");
}
