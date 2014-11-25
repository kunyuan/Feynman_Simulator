//
//  weight_MonteCarlo.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/23/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "weight_MonteCarlo.h"
using namespace mc::weight;

G::G(model Model, const Lattice &lat, real beta,
     const std::vector<real> &hopping,
     const std::vector<real> &RealChemicalPotential,
     real ExternalField, weight0::TauSymmetry TauSymmetry)
    : weight0::G(Model, lat, beta, hopping, RealChemicalPotential, ExternalField, TauSymmetry)
{
    _MeasureWeight.Allocate(GetShape());
    //initialize _MeasureWeight to an unit function
    _MeasureWeight = Complex(1.0, 0.0);
}

Complex G::Weight(const Site &rin, const Site &rout, real tin, real tout, spin SpinIn, spin SpinOut, bool IsMeasure)
{
    auto dist = _Lat.Dist(rin, rout);
    if (IsMeasure)
        return _MeasureWeight[_Map.SpinIndex(SpinIn, SpinOut)]
                             [_Map.SublatIndex(dist)]
                             [_Map.CoordiIndex(dist)]
                             [_Map.TauIndex(tin, tout)];
    else
        return _TauSymmetryFactor *
               _SmoothTWeight[_Map.SpinIndex(SpinIn, SpinOut)]
                             [_Map.SublatIndex(dist)]
                             [_Map.CoordiIndex(dist)]
                             [_Map.TauIndex(tin, tout)];
}

Complex G::Weight(int dir, const Site &r1, const Site &r2, real t1, real t2, spin Spin1, spin Spin2, bool IsMeasure)
{
    if (dir == IN) {
        auto dist = _Lat.Dist(r1, r2);

        if (IsMeasure)
            return _MeasureWeight[_Map.SpinIndex(Spin1, Spin2)]
                                 [_Map.SublatIndex(dist)]
                                 [_Map.CoordiIndex(dist)]
                                 [_Map.TauIndex(t1, t2)];

        else
            return _TauSymmetryFactor *
                   _SmoothTWeight[_Map.SpinIndex(Spin1, Spin2)]
                                 [_Map.SublatIndex(dist)]
                                 [_Map.CoordiIndex(dist)]
                                 [_Map.TauIndex(t1, t2)];
    }
    else {
        auto dist = _Lat.Dist(r2, r1);
        if (IsMeasure)
            return _MeasureWeight[_Map.SpinIndex(Spin2, Spin1)]
                                 [_Map.SublatIndex(dist)]
                                 [_Map.CoordiIndex(dist)]
                                 [_Map.TauIndex(t2, t1)];

        else
            return _TauSymmetryFactor *
                   _SmoothTWeight[_Map.SpinIndex(Spin2, Spin1)]
                                 [_Map.SublatIndex(dist)]
                                 [_Map.CoordiIndex(dist)]
                                 [_Map.TauIndex(t2, t1)];
    }
}
