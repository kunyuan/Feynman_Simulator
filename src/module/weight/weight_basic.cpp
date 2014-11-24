//
//  weight_basic.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/21/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "weight_basic.h"
#include <iostream>
#include "utility/logger.h"
#include <math.h>

using namespace std;
using namespace weight0;

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

Basic::Basic(const Lattice &lat, real beta, SpinNum spin_num, model Model,
             TauSymmetry Symmetry, string name)
    : _Lat(lat),
      _Beta(beta),
      _TauSymmetryFactor(int(Symmetry)),
      _Name(name),
      _Model(Model),
      _SpinNum(int(spin_num))

{
    _dBeta = beta / MAX_TAU_BIN;
    _dBetaInverse = 1.0 / _dBeta;
    CheckVec2Index(lat);
}

vector<uint> Basic::GetShape()
{
    auto SpinVol = static_cast<uint>(pow(2, _SpinNum));
    vector<uint> shape({SpinVol, (uint)_Lat.SublatVol2, (uint)_Lat.Vol});
    shape.push_back(MAX_TAU_BIN);
    return shape;
}

vector<uint> Basic::GetSpaceShape()
{
    vector<uint> shape;
    for (auto e : _Lat.Size)
        shape.push_back(e);
    return shape;
}

vector<uint> Basic::GetSpaceTimeShape()
{
    vector<uint> shape = GetSpaceShape();
    shape.push_back(MAX_TAU_BIN);
    return shape;
}

void Basic::Reset(real beta)
{
    _Beta = beta;
    _dBeta = beta / MAX_TAU_BIN;
    _dBetaInverse = 1.0 / _dBeta;
    //TODO: please implement how to reset the weight here
}

int Basic::SublatIndex(const Distance &dist)
{
    return dist.SublatIndex;
}

int Basic::CoordiIndex(const Distance &dist)
{
    return dist.CoordiIndex;
}

int Basic::TauIndex(real tau)
{
    if (DEBUGMODE && tau < -_Beta || tau >= _Beta)
        LOG_INFO("tau=" << tau << " is out of the range ["
                        << -_Beta << "," << _Beta << ")");
    //TODO: mapping between tau and bin

    int bin = tau < 0 ? floor(tau * _dBetaInverse) + MAX_TAU_BIN
                      : floor(tau * _dBetaInverse);
    if (DEBUGMODE && bin < 0 || tau >= MAX_TAU_BIN) {
        LOG_INFO("tau=" << tau << " is out of the range ["
                        << -_Beta << "," << _Beta << ")");
        LOG_INFO("bin=" << bin << " is out of the range ["
                        << 0 << "," << MAX_TAU_BIN << "]");
    }
    return bin;
}

int Basic::TauIndex(real t_in, real t_out)
{
    return TauIndex(t_out - t_in);
}

real Basic::IndexToTau(int Bin)
{
    //TODO: mapping between tau and bin
    return Bin * _dBeta + _Beta / 2;
}

BasicWithTwoSpins::BasicWithTwoSpins(const Lattice &lat, real Beta, model Model,
                                     TauSymmetry Symmetry, std::string Name)
    : Basic(lat, Beta, TwoSpins, Model, Symmetry, Name)
{
}

int BasicWithTwoSpins::SpinIndex(spin SpinIn, spin SpinOut)
{
    return SpinIn * SPIN + SpinOut;
}

bool BasicWithTwoSpins::IsSameSpin(int spindex)
{
    return (spindex == 0 || spindex == 2);
}

BasicWithFourSpins::BasicWithFourSpins(const Lattice &lat, real Beta, model Model,
                                       TauSymmetry Symmetry, std::string Name)
    : Basic(lat, Beta, FourSpins, Model, Symmetry, Name)
{
}
//First In/Out: direction of WLine; Second In/Out: direction of Vertex
int BasicWithFourSpins::SpinIndex(spin SpinInIn, spin SpinInOut, spin SpinOutIn, spin SpinOutOut)
{
    return SpinInIn * SPIN3 + SpinInOut * SPIN2 +
           SpinOutIn * SPIN + SpinOutOut;
}
int BasicWithFourSpins::SpinIndex(spin *TwoSpinIn, spin *TwoSpinOut)
{
    return SpinIndex(TwoSpinIn[0], TwoSpinIn[1],
                     TwoSpinOut[0], TwoSpinOut[1]);
}

std::vector<int> BasicWithFourSpins::GetSpinIndexVector(SpinFilter filter)
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
