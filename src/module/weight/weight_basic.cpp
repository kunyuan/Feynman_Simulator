//
//  weight_basic.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/21/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "weight_basic.h"
#include "utility/logger.h"
#include "utility/cnpy.h"
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

Basic::Basic(const Lattice &lat, real beta, SpinNum spin_num,
             TauSymmetry Symmetry, string name)
    : _Lat(lat),
      _Beta(beta),
      _TauSymmetryFactor(int(Symmetry)),
      _Name(name),
      _SpinNum(int(spin_num))
{
    _dBeta = beta / MAX_TAU_BIN;
    _dBetaInverse = 1.0 / _dBeta;
    CheckVec2Index(lat);

    auto SpinVol = static_cast<uint>(pow(2, _SpinNum));
    _Shape = vector<uint>({SpinVol, (uint)_Lat.SublatVol2, (uint)_Lat.Vol, MAX_TAU_BIN});
    for (auto e : _Lat.Size)
        _SpaceTimeShape.push_back(e);
    _SpaceTimeShape.push_back(MAX_TAU_BIN);

    _SmoothTWeight.Allocate(GetShape());
    _DeltaTWeight.Allocate(GetShape());
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
    _dBeta = beta / MAX_TAU_BIN;
    _dBetaInverse = 1.0 / _dBeta;
    //TODO: please implement how to reset the weight here
}

int Basic::GetTauSymmetryFactor(real t_in, real t_out)
{
    return _TauSymmetryFactor * ((t_out > t_in) ? 1 : -1);
}

void Basic::BuildNew(model Model)
{
    _Model = Model;
}

void Basic::BuildTest()
{
    _Model = model::TEST;
}

const string SMOOTH = ".SmoothT";
const string DELTA = ".DeltaT";

template <typename T>
bool LoadMatrix(T &matrix, const string &FileName, const string &Name)
{
    cnpy::NpyArray weight = cnpy::npz_load(FileName, Name);
    if (weight.data == nullptr) {
        ABORT("Can't find " << Name << ".Weight in .npz data file!");
        return false;
    }
    //assignment here will copy data in weight.data into *this
    matrix = (reinterpret_cast<Complex *>(weight.data));
    return true;
}

template <typename T>
void SaveMatrix(T &matrix, const string &FileName, const std::string Mode,
                const string &Name, uint *Shape, int Dim)
{
    cnpy::npz_save(FileName, Name, matrix(), Shape, Dim, Mode);
}

bool Basic::Load(const std::string &FileName)
{
    return LoadMatrix(_SmoothTWeight, FileName, _Name + SMOOTH) &&
           LoadMatrix(_DeltaTWeight, FileName, _Name + DELTA);
}
void Basic::Save(const std::string &FileName, const std::string Mode)
{
    SaveMatrix(_SmoothTWeight, FileName, Mode, _Name + SMOOTH, GetShape(), 4);
    SaveMatrix(_DeltaTWeight, FileName, "a", _Name + DELTA, GetShape(), 3);
}
