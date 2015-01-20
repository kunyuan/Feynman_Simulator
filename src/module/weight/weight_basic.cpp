//
//  weight_basic.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/21/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "weight_basic.h"
#include "utility/logger.h"
#include "utility/dictionary.h"
#include <math.h>

using namespace std;
using namespace weight;

Basic::Basic(const Lattice& lat, real beta, uint MaxTauBin, SpinNum spin_num,
             TauSymmetry Symmetry, string name)
    : _Lat(lat)
    , _Beta(beta)
    , _TauSymmetryFactor(int(Symmetry))
    , _Name(name)
    , _SpinNum(int(spin_num))
{
    _MaxTauBin = MaxTauBin;
    _dBeta = beta / _MaxTauBin;
    _dBetaInverse = 1.0 / _dBeta;

    auto SpinVol = static_cast<uint>(pow(2, _SpinNum));
    _Shape = vector<uint>({ SpinVol, (uint)_Lat.SublatVol2, (uint)_Lat.Vol, _MaxTauBin });
}

uint* Basic::GetShape()
{
    return _Shape.data();
}

void Basic::Reset(real beta)
{
    _Beta = beta;
    _dBeta = beta / _MaxTauBin;
    _dBetaInverse = 1.0 / _dBeta;
    //TODO: please implement how to reset the weight here
}

int Basic::GetTauSymmetryFactor(real t_in, real t_out) const
{
    return (t_out > t_in) ? 1 : _TauSymmetryFactor;
}

const string SMOOTH = "SmoothT";
const string DELTA = "DeltaT";

bool DeltaTArray::FromDict(const Dictionary& dict)
{
    auto arr = dict.Get<Python::ArrayObject>(DELTA);
    ASSERT_ALLWAYS(Equal(arr.Shape().data(), GetShape(), 3), "Shape should match!");
    Assign(arr.Data<Complex>());
    return true;
}
Dictionary DeltaTArray::ToDict()
{
    Dictionary dict;
    dict[DELTA] = Python::ArrayObject(Data(), GetShape(), 3);
    return dict;
}
bool SmoothTArray::FromDict(const Dictionary& dict)
{
    auto arr = dict.Get<Python::ArrayObject>(SMOOTH);
    ASSERT_ALLWAYS(Equal(arr.Shape().data(), GetShape(), 4), "Shape should match!");
    Assign(arr.Data<Complex>());
    return true;
}
Dictionary SmoothTArray::ToDict()
{
    Dictionary dict;
    dict[SMOOTH] = Python::ArrayObject(Data(), GetShape(), 4);
    return dict;
}
