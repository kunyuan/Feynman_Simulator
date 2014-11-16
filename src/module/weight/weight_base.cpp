//
//  observable.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/13/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "weight_base.h"
#include "utility/cnpy.h"
#include "utility/abort.h"
#include "utility/scopeguard.h"

using namespace std;
using namespace Array;
using namespace weight;

WeightNoMeasure::WeightNoMeasure(const Lattice &lat, real beta,
                                 int order, bool IsTauSymmetric_,
                                 int SpinVol, string name)
    : _Lat(lat),
      _Beta(beta),
      _Order(order),
      _IsTauSymmetric(IsTauSymmetric_),
      _Name(name)

{
    _dBeta = beta / MAX_BIN;
    _dBetaInverse = 1.0 / _dBeta;

    _Shape[ORDER] = order;
    _Shape[SP] = SpinVol;
    _Shape[SUB] = lat.SublatVol * lat.SublatVol;
    _Shape[VOL] = lat.Vol;
    _Shape[TAU] = MAX_BIN;
    //if you want to change the order of _Shape, don't forget to take care of Dyson module

    SmoothWeight.Allocate(Shape());
    SmoothWeight = 0.0;
    //use _Shape[SP] to _Shape[TAU] to construct array4

    DeltaTWeight.Allocate(Shape());
    DeltaTWeight = 0.0;
    //use _Shape[SP] to _Shape[Vol] to construct array3

    for (int i = 0; i < lat.Dimension; i++)
        _SpaceTimeShape[i] = lat.Size[i];
    _SpaceTimeShape[lat.Dimension] = MAX_BIN;
    //Lx*Ly*Lz*Lt

    if (!_CheckVec2Index()) {
        string message = "The mapping between vector and index of \
        the lattice should look like\n (i,j,k) -> k + n3*j + n2*n3*i; \
        \n n1, n2, n3 : dimensions in three directions; \
        which is required by fft on spatial dimensions";
        ABORT(message);
    }
}

unsigned int *WeightNoMeasure::Shape()
{
    return _Shape + SP;
}

void WeightNoMeasure::Reset(real beta)
{
    _Beta = beta;
    _dBeta = beta / MAX_BIN;
    _dBetaInverse = 1.0 / _dBeta;
    //TODO: please implement how to reset the weight here
}

int WeightNoMeasure::SpinIndex(spin SpinIn, spin SpinOut)
{
    return SpinIn * SPIN + SpinOut;
}

int WeightNoMeasure::SpinIndex(spin SpinInIn, spin SpinInOut, spin SpinOutIn, spin SpinOutOut)
{
    return SpinInIn * SPIN3 + SpinInOut * SPIN2 +
           SpinOutIn * SPIN + SpinOutOut;
}

int WeightNoMeasure::SpinIndex(spin *TwoSpinIn, spin *TwoSpinOut)
{
    return TwoSpinIn[0] * SPIN3 + TwoSpinIn[1] * SPIN2 +
           TwoSpinOut[0] * SPIN + TwoSpinOut[1];
}

vector<int> WeightNoMeasure::GetSpinIndexVector_FourSpinsFileter(SpinFilter filter)
{
    vector<int> list;
    for (int InIn = 0; InIn < 2; InIn++)
        for (int InOut = 0; InOut < 2; InOut++)
            for (int OutIn = 0; OutIn < 2; OutIn++)
                for (int OutOut = 0; OutOut < 2; OutOut++) {
                    bool flag = false;
                    if (filter == UpUp2UpUp &&
                        InIn == InOut && InIn == OutIn && InIn == OutOut)
                        flag = true;
                    if (filter == UpDown2UpDown &&
                        InIn == InOut && OutIn == OutOut && InIn == FLIP(OutIn))
                        flag = true;
                    if (filter == UpDown2DownUp &&
                        InIn == FLIP(InOut) && OutIn == FLIP(OutOut) && InIn == FLIP(OutIn))
                        flag = true;
                    if (flag)
                        list.push_back(SpinIndex(spin(InIn), spin(InOut),
                                                 spin(OutIn), spin(OutOut)));
                }
    return list;
}

int WeightNoMeasure::TauToBin(real tau)
{
    if (DEBUGMODE && tau < -_Beta || tau >= _Beta)
        LOG_INFO("tau=" << tau << " is out of the range ["
                        << -_Beta << "," << _Beta << ")");
    //TODO: mapping between tau and bin

    int bin = tau < 0 ? floor(tau * _dBetaInverse) + MAX_BIN
                      : floor(tau * _dBetaInverse);
    if (DEBUGMODE && bin < 0 || tau >= MAX_BIN) {
        LOG_INFO("tau=" << tau << " is out of the range ["
                        << -_Beta << "," << _Beta << ")");
        LOG_INFO("bin=" << bin << " is out of the range ["
                        << 0 << "," << MAX_BIN << "]");
    }
    return bin;
}

int WeightNoMeasure::TauSymmetry(real t_in, real t_out)
{
    return (_IsTauSymmetric || t_out > t_in) ? 1 : -1;
}

int WeightNoMeasure::TauToBin(real t_in, real t_out)
{
    return TauToBin(t_out - t_in);
}

real WeightNoMeasure::BinToTau(int Bin)
{
    //TODO: mapping between tau and bin
    return Bin * _dBeta + _Beta / 2;
}

void WeightNoMeasure::Save(const std::string &FileName, std::string Mode)
{
    cnpy::npz_save(FileName, _Name + ".Weight", SmoothWeight(), Shape(), 4, Mode);
    cnpy::npz_save(FileName, _Name + ".DeltaTWeight", DeltaTWeight(), Shape(), 3, "a");
}

bool WeightNoMeasure::Load(const std::string &FileName)
{
    cnpy::npz_t NpzMap = cnpy::npz_load(FileName);
    ON_SCOPE_EXIT([&] {NpzMap.destruct(); });

    cnpy::NpyArray weight = NpzMap[_Name + ".Weight"];
    if (weight.data == nullptr)
        ABORT("Can't find " << _Name << ".Weight in .npz data file!");
    //assignment here will copy data in weight.data into _Weight
    SmoothWeight = reinterpret_cast<Complex *>(weight.data);

    cnpy::NpyArray delta_weight = NpzMap[_Name + ".DeltaTWeight"];
    if (delta_weight.data == nullptr)
        ABORT("Can't find " << _Name << ".DeltaTWeight in .npz data file!");
    //assignment here will copy data in weight.data into _Weight
    DeltaTWeight = reinterpret_cast<Complex *>(delta_weight.data);
    return true;
}
/**
*  Check if the mapping between vector and index of the lattice looks like
*            (i,j,k) -> k + n3*j + n2*n3*i;
*     n1, n2, n3 : dimensions in three directions;
*   which is required by fft on spatial dimensions
*/
bool WeightNoMeasure::_CheckVec2Index()
{
    Vec<int> v;
    for (int index = 0; index < _Lat.Vol; index++) {
        int j = index;
        for (int i = D - 1; i > 0; i--) {
            v[i] = j % _Lat.Size[i];
            j /= _Lat.Size[i];
        }
        v[0] = j;
        if (v != _Lat.Index2Vec(index))
            return false;
    }
    return true;
}
