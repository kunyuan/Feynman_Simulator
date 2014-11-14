//
//  weight_IO.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/20/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "weight_inherit.h"

using namespace std;
using namespace Array;
using namespace weight;

Sigma::Sigma(const Lattice &lat, real beta, int order, bool IsTauSymmetric)
    : WeightNeedMeasure(lat, beta, order,
                        IsTauSymmetric, SPIN2, "Sigma", Norm::Weight())
{
}

Complex Sigma::Weight(const Site &rin, const Site &rout, real tin, real tout, spin SpinIn, spin SpinOut)
{
    auto dist = _Lat.Dist(rin, rout);
    return TauSymmetry(tin, tout) *
           SmoothWeight[SpinIndex(SpinIn, SpinOut)]
                       [dist.SublatIndex]
                       [dist.CoordiIndex]
                       [TauToBin(tin, tout)];
}

Complex Sigma::WeightOfDelta(spin SpinIn, spin SpinOut)
{
    //TODO: add Delta expression of Sigma here, you may need to add new API to measure it
    return 0.0;
}

void Sigma::Measure(const Site &rin, const Site &rout, real tin, real tout, spin SpinIn, spin SpinOut, int order, const Complex &weight)
{
    if (DEBUGMODE && order <= 0)
        LOG_ERROR("Too small order=" << order);
    int spin = SpinIndex(SpinIn, SpinOut);
    auto dist = _Lat.Dist(rin, rout);
    int tau = TauToBin(tin, tout);
    _WeightAccu[order - 1][spin][dist.SublatIndex]
               [dist.CoordiIndex][tau] += weight * TauSymmetry(tin, tout);
    if (spin == 0 && dist.SublatIndex == 0 && dist.CoordiIndex == 0 && tau == 0)
        _WeightErrorEstimator[order - 1].Measure(weight * TauSymmetry(tin, tout));
}

/************************   Polarization   *********************************/
//
Polar::Polar(const Lattice &lat, real beta, int order)
    : WeightNeedMeasure(lat, beta, order, false, SPIN4, "Polar", Norm::Weight())
{
}

Complex Polar::Weight(const Site &rin, const Site &rout, real tin, real tout, spin *SpinIn, spin *SpinOut)
{
    auto dist = _Lat.Dist(rin, rout);
    return SmoothWeight[SpinIndex(SpinIn, SpinOut)]
                       [dist.SublatIndex]
                       [dist.CoordiIndex]
                       [TauToBin(tout - tin)];
}

void Polar::Measure(const Site &rin, const Site &rout, real tin, real tout, spin *SpinIn, spin *SpinOut, int order, const Complex &weight)
{
    if (DEBUGMODE && order <= 0)
        LOG_ERROR("Too small order=" << order);
    int spin = SpinIndex(SpinIn, SpinOut);
    auto dist = _Lat.Dist(rin, rout);
    int tau = TauToBin(tin, tout);
    _WeightAccu[order - 1][spin][dist.SublatIndex][dist.CoordiIndex][tau] += weight;
    if (spin == 0 && dist.SublatIndex == 0 && dist.CoordiIndex == 0 && tau == 0)
        _WeightErrorEstimator[order - 1].Measure(weight);
}

/***********************  G  **************************************/

G::G(const Lattice &lat, real beta, int order, real ExternalField, bool IsTauSymmetric)
    : WeightNoMeasure(lat, beta, order, IsTauSymmetric, SPIN4, "G")
{
    _ExternalField = ExternalField;
    BareWeight.Allocate(Shape());
    //use _Shape[SP] to _Shape[TAU] to construct array3
    InitialWithBare();
}

Complex G::Weight(const Site &rin, const Site &rout, real tin, real tout, spin SpinIn, spin SpinOut, bool IsMeasure)
{
    auto dist = _Lat.Dist(rin, rout);
    if (IsMeasure)
        //TODO: define the measuring weight for G
        return Complex(1.0, 0.0);
    else
        return TauSymmetry(tin, tout) *
               SmoothWeight[SpinIndex(SpinIn, SpinOut)]
                           [dist.SublatIndex]
                           [dist.CoordiIndex]
                           [TauToBin(tin, tout)];
}

Complex G::Weight(int dir, const Site &r1, const Site &r2, real t1, real t2, spin Spin1, spin Spin2, bool IsMeasure)
{
    if (IsMeasure)
        //TODO: define the measuring weight for G
        return Complex(1.0, 0.0);

    else if (dir == IN) {
        auto dist = _Lat.Dist(r1, r2);
        return TauSymmetry(t1, t2) *
               SmoothWeight[SpinIndex(Spin1, Spin2)]
                           [dist.SublatIndex]
                           [dist.CoordiIndex]
                           [TauToBin(t1, t2)];
    }
    else {
        auto dist = _Lat.Dist(r2, r1);
        return TauSymmetry(t2, t1) *
               SmoothWeight[SpinIndex(Spin2, Spin1)]
                           [dist.SublatIndex]
                           [dist.CoordiIndex]
                           [TauToBin(t2, t1)];
    }
}

/***********************  W  **************************************/

W::W(const Lattice &lat, real Beta, int order,
     const vector<real> &Interaction_, real ExternalField)
    : WeightNoMeasure(lat, Beta, order, false, SPIN4, "W")
{
    _Interaction = Interaction_;
    _ExternalField = ExternalField;
    BareWeight.Allocate(Shape());
    //use _Shape[SP] to _Shape[VOL] to construct array3
    InitialWithBare();
}

Complex W::Weight(const Site &rin, const Site &rout, real tin, real tout, spin *SpinIn, spin *SpinOut, bool IsWorm, bool IsMeasure, bool IsDelta)
{
    if (IsMeasure)
        //TODO: define the measuring weight of W
        return Complex(1.0, 0.0);

    auto distance = _Lat.Dist(rin, rout);
    if (IsDelta)
        //TODO: define the delta function here! IsWorm==true and IsWorm==false
        return DeltaTWeight[SpinIndex(SpinIn, SpinOut)]
                           [distance.SublatIndex]
                           [distance.CoordiIndex];

    else if (IsWorm)
        //define your fake function here
        return SmoothWeight[SpinIndex(UP, UP)]
                           [distance.SublatIndex]
                           [distance.CoordiIndex]
                           [TauToBin(tin, tout)];
    else
        return SmoothWeight[SpinIndex(SpinIn, SpinOut)]
                           [distance.SublatIndex]
                           [distance.CoordiIndex]
                           [TauToBin(tin, tout)];
}

Complex W::Weight(int dir, const Site &r1, const Site &r2, real t1, real t2, spin *Spin1, spin *Spin2, bool IsWorm, bool IsMeasure, bool IsDelta)
{
    int spinindex, subindex, coordindex, tau;
    if (dir == IN) {
        auto dist = _Lat.Dist(r1, r2);
        spinindex = SpinIndex(Spin1, Spin2);
        subindex = dist.SublatIndex;
        coordindex = dist.CoordiIndex;
        tau = TauToBin(t1, t2);
    }
    else {
        auto dist = _Lat.Dist(r2, r1);
        spinindex = SpinIndex(Spin2, Spin1);
        subindex = dist.SublatIndex;
        coordindex = dist.CoordiIndex;
        tau = TauToBin(t2, t1);
    }

    if (IsMeasure)
        //TODO: define the measuring weight of W
        return Complex(1.0, 0.0);

    if (IsDelta)
        //TODO: define the delta function here! IsWorm==true and IsWorm==false
        return DeltaTWeight[spinindex]
                           [subindex]
                           [coordindex];
    else if (IsWorm)
        //define your fake function here
        return SmoothWeight[SpinIndex(UP, UP)][subindex][coordindex][tau];
    else
        return SmoothWeight[spinindex][subindex][coordindex][tau];
}

void W::WriteBareToASCII()
{
    _Lat.PlotLattice();
    Vec<int> offset;
    for (int i = 0; i < D; i++)
        offset[i] = _Lat.Size[i] / 2 - 1;

    ofstream os("interaction.py", ios::out);
    int spin_index = SpinIndex(UP, UP);
    os << "line=[" << endl;
    for (int sub = 0; sub < _Shape[SUB]; sub++)
        for (int coord = 0; coord < _Shape[VOL]; coord++) {
            real bare_weight = mod(BareWeight[spin_index][sub][coord]);
            if (!Zero(bare_weight)) {
                Distance dis(sub, coord);
                Site start = _Lat.GetSite(dis, IN);
                Site end = _Lat.GetSite(dis, OUT);
                os << "[" << _Lat.GetRealVec(start, offset).PrettyString() << ","
                   << _Lat.GetRealVec(end, offset).PrettyString() << ","
                   << start.Sublattice << "]," << endl;
            }
        }
    os << "]" << endl;
}
