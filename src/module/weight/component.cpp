//
//  component.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/22/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "component.h"
#include "weight_initializer.h"
#include "utility/cnpy.h"
#include <tuple>

using namespace weight0;
using namespace std;

G::G(model Model, const Lattice &lat, real beta,
     const std::vector<real> &Hopping,
     const std::vector<real> &RealChemicalPotential,
     real ExternalField, TauSymmetry Symmetry)
    : weight0::Basic(Model, lat, beta, SPIN2, Symmetry, "G"),
      _Map(IndexMapSPIN2(beta, lat))
{
    _Hopping = Hopping;
    _ExternalField = ExternalField;
    _RealChemicalPotential = RealChemicalPotential;
    //use _Shape[SP] to _Shape[TAU] to construct array3
    _MeasureWeight.Allocate(GetShape());
    //initialize _MeasureWeight to an unit function
    _MeasureWeight = Complex(1.0, 0.0);
}

void G::BuildNew()
{
    GInitializer(*this).BuildNew();
}

void G::BuildTest()
{
    GInitializer(*this).BuildTest();
}

void G::Reset(real Beta)
{
    Basic::Reset(Beta);
    _Map = IndexMapSPIN2(Beta, _Lat);
}

W::W(model Model, const Lattice &lat, real Beta,
     const vector<real> &Interaction_, real ExternalField)
    : weight0::Basic(Model, lat, Beta, SPIN4, TauSymmetric, "W"),
      _Map(IndexMapSPIN4(Beta, lat))
{
    _Interaction = Interaction_;
    _ExternalField = ExternalField;
    //use _Shape[SP] to _Shape[VOL] to construct array3
    _MeasureWeight.Allocate(GetShape());
    //initialize _MeasureWeight to an unit function
    _MeasureWeight = Complex(1.0, 0.0);
}

void W::BuildNew()
{
    WInitializer(*this).BuildNew();
}

void W::BuildTest()
{
    WInitializer(*this).BuildTest();
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

Sigma::Sigma(model Model, const Lattice &lat, real Beta, int MaxOrder, TauSymmetry Symmetry)
    : weight0::Basic(Model, lat, Beta, SPIN2, Symmetry, "Sigma"),
      _Map(IndexMapSPIN2(Beta, lat)),
      Estimator(Beta, MaxOrder, "Sigma", Norm::Weight(), GetShape())
{
}

void Sigma::Reset(real Beta)
{
    Basic::Reset(Beta);
    _Map = IndexMapSPIN2(Beta, _Lat);
    Estimator.ReWeight(Beta);
}

bool Sigma::Load(const std::string &FileName)
{
    return Estimator.Load(FileName) &&
           Basic::Load(FileName);
}
void Sigma::Save(const std::string &FileName, const std::string Mode)
{
    Estimator.Save(FileName, Mode);
    Basic::Save(FileName, "a");
}

Polar::Polar(model Model, const Lattice &lat, real Beta, int MaxOrder)
    : weight0::Basic(Model, lat, Beta, SPIN4, TauSymmetric, "Polar"),
      _Map(IndexMapSPIN4(Beta, lat)),
      Estimator(Beta, MaxOrder, "Polar", Norm::Weight(), GetShape())
{
}

void Polar::Reset(real Beta)
{
    Basic::Reset(Beta);
    _Map = IndexMapSPIN4(Beta, _Lat);
    Estimator.ReWeight(Beta);
}

bool Polar::Load(const std::string &FileName)
{
    return Estimator.Load(FileName) &&
           Basic::Load(FileName);
}
void Polar::Save(const std::string &FileName, const std::string Mode)
{
    Estimator.Save(FileName, Mode);
    Basic::Save(FileName, "a");
}
