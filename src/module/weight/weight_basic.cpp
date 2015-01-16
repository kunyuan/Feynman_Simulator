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

/**
*  Check if the mapping between vector and index of the lattice looks like
*            (i,j,k) -> k + n3*j + n2*n3*i;
*     n1, n2, n3 : dimensions in three directions;
*   which is required by fft on spatial dimensions
*/
void CheckVec2Index(Lattice _Lat)
{
    Vec<int> v;
    bool flag = true;
    for (int index = 0; index < _Lat.Vol; index++) {
        int j = index;
        for (int i = D - 1; i > 0; i--) {
            v[i] = j % _Lat.Size[i];
            j /= _Lat.Size[i];
        }
        v[0] = j;
        if (v != _Lat.Index2Vec(index)) {
            flag = false;
            break;
        }
    }
    if (!flag) {
        string message = "The mapping between vector and index of \
        the lattice should look like\n (i,j,k) -> k + n3*j + n2*n3*i; \
        \n n1, n2, n3 : dimensions in three directions; \
        which is required by fft on spatial dimensions";
        ABORT(message);
    }
}

Basic::Basic(const Lattice &lat, real beta, uint MaxTauBin, SpinNum spin_num,
             TauSymmetry Symmetry, string name)
    : _Lat(lat), _Beta(beta), _TauSymmetryFactor(int(Symmetry)), _Name(name), _SpinNum(int(spin_num))
{
    _MaxTauBin = MaxTauBin;
    _dBeta = beta / _MaxTauBin;
    _dBetaInverse = 1.0 / _dBeta;
    CheckVec2Index(lat);

    auto SpinVol = static_cast<uint>(pow(2, _SpinNum));
    _Shape = vector<uint>({SpinVol, (uint)_Lat.SublatVol2, (uint)_Lat.Vol, _MaxTauBin});
    for (auto e : _Lat.Size)
        _SpaceTimeShape.push_back(e);
    _SpaceTimeShape.push_back(_MaxTauBin);

    _SmoothTWeight.Allocate(GetShape());
    _SmoothTWeight = Complex(0.0, 0.0);
    _DeltaTWeight.Allocate(GetShape());
    _DeltaTWeight = Complex(0.0, 0.0);
}

uint *Basic::GetShape()
{
    return _Shape.data();
}

uint *Basic::GetSpaceTimeShape()
{
    return _SpaceTimeShape.data();
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

bool Basic::FromDict(const Dictionary &dict)
{
    bool flagSmooth = dict.HasKey(SMOOTH);
    bool flagDelta = dict.HasKey(DELTA);
    if (flagSmooth) {
        auto arr = dict.Get<Python::ArrayObject>(SMOOTH);
        ASSERT_ALLWAYS(Equal(arr.Shape().data(), GetShape(), 4), "Shape should match!");
        _SmoothTWeight = arr.Data<Complex>();
    }
    if (flagDelta) {
        auto arr = dict.Get<Python::ArrayObject>(DELTA);
        ASSERT_ALLWAYS(Equal(arr.Shape().data(), GetShape(), 3), "Shape should match!");
        _DeltaTWeight = arr.Data<Complex>();
    }
    ASSERT_ALLWAYS(flagSmooth || flagDelta, "Come on! Neither SmoothT nor DeltaT array exist!");
    return true;
}

Dictionary Basic::ToDict()
{
    Dictionary dict;
    dict[SMOOTH] = Python::ArrayObject(_SmoothTWeight(), GetShape(), 4);
    dict[DELTA] = Python::ArrayObject(_DeltaTWeight(), GetShape(), 3);
    return dict;
}
