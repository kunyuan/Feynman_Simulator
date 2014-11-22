//
//  weight_array.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/20/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "weight_matrix.h"
#include "utility/cnpy.h"
using namespace std;
using namespace weight;

vector<uint> GetShape(Lattice lat, SpinNum spin_num, bool HasTime)
{
    auto SpinVol = static_cast<uint>(pow(2, spin_num));
    vector<uint> shape({SpinVol, (uint)lat.SublatVol2, (uint)lat.Vol});
    if (HasTime)
        shape.push_back(MAX_TAU_BIN);
    return shape;
}

vector<uint> GetSpaceShape(Lattice lat)
{
    vector<uint> shape;
    for (auto e : lat.Size)
        shape.push_back(e);
    return shape;
}

vector<uint> GetSpaceTimeShape(Lattice lat)
{
    vector<uint> shape = GetSpaceShape(lat);
    shape.push_back(MAX_TAU_BIN);
    return shape;
}

/*********************************************************************/
/*                          SmoothTWeight                             */
/*********************************************************************/
SmoothTWeight::SmoothTWeight(real Beta, Lattice lat, SpinNum spin_num, bool IsTauSymmetric,
                             const std::string Name)
{
    _Name = Name;
    _Beta = Beta;
    _IsTauSymmetric = IsTauSymmetric;
    Shape = GetShape(lat, spin_num, true);
    SpaceTimeShape = GetSpaceTimeShape(lat);
}

void SmoothTWeight::Activate()
{
    Allocate(Shape.data());
}

bool SmoothTWeight::Load(const std::string &FileName)
{
    cnpy::NpyArray weight = cnpy::npz_load(FileName, _Name + ".Weight");
    if (weight.data == nullptr) {
        ABORT("Can't find " << _Name << ".Weight in .npz data file!");
        return false;
    }
    //assignment here will copy data in weight.data into *this
    Array::array4<Complex>::operator=(reinterpret_cast<Complex *>(weight.data));
    return true;
}

void SmoothTWeight::Save(const std::string &FileName, const std::string Mode)
{
    cnpy::npz_save(FileName, _Name + ".Weight", (*this)(), Shape.data(), 4, Mode);
}

/**
*  fft::FORTH would change antisymmetric _Weight to be symmetric, and always performed in Tau space
*  fft::BACK would change symmetric _Weight to be antisymmetric, and always performed in Omega space
*/
void SmoothTWeight::_ChangeSymmetry(fft::Dir direction)
{
    int sign = static_cast<int>(direction);
    Complex PhaseFactor[MAX_TAU_BIN];
    for (int i = 0; i < Shape[TAU]; i++)
        PhaseFactor[i] = exp(Complex(0.0, sign * BinToTau(i) * PI));

    int NumOfTimeSeries = Shape[SP] * Shape[SUB] * Shape[VOL];
    for (int i = 0; i < NumOfTimeSeries; i += Shape[TAU])
        for (int j = 0; j < Shape[TAU]; j++)
            SmoothWeight(i * Shape[TAU] + j) *= PhaseFactor[j];
}

/*********************************************************************/
/*                          DeltaTWeight                             */
/*********************************************************************/
DeltaTWeight::DeltaTWeight(Lattice lat, SpinNum spin_num,
                           const std::string Name)
{
    _Name = Name;
    Shape = GetShape(lat, spin_num, false);
    SpaceShape = GetSpaceShape(lat);
}

void DeltaTWeight::Activate()
{
    Allocate(Shape.data());
}

bool DeltaTWeight::Load(const std::string &FileName)
{
    cnpy::NpyArray weight = cnpy::npz_load(FileName, _Name + ".Weight");
    if (weight.data == nullptr) {
        ABORT("Can't find " << _Name << ".Weight in .npz data file!");
        return false;
    }
    //assignment here will copy data in weight.data into *this
    Array::array3<Complex>::operator=(reinterpret_cast<Complex *>(weight.data));
    return true;
}

void DeltaTWeight::Save(const std::string &FileName, const std::string Mode)
{
    cnpy::npz_save(FileName, _Name + ".Weight", (*this)(), Shape.data(), 3, Mode);
}
