//
//  weight_basic.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/21/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "weight_basic.h"
#include "weight_matrix.h"
#include <iostream>
#include "utility/logger.h"
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

Basic::Basic(const Lattice &lat, real beta, model Model,
             bool IsTauSymmetric_, string name)
    : _Lat(lat),
      _Beta(beta),
      _IsTauSymmetric(IsTauSymmetric_),
      _Name(name),
      _Model(Model)

{
    _dBeta = beta / MAX_TAU_BIN;
    _dBetaInverse = 1.0 / _dBeta;
    CheckVec2Index(lat);
}

int Basic::SpinIndex(spin SpinIn, spin SpinOut)
{
    return SpinIn * SPIN + SpinOut;
}

int Basic::SpinIndex(spin SpinInIn, spin SpinInOut, spin SpinOutIn, spin SpinOutOut)
{
    return SpinInIn * SPIN3 + SpinInOut * SPIN2 +
           SpinOutIn * SPIN + SpinOutOut;
}

int Basic::SpinIndex(spin *TwoSpinIn, spin *TwoSpinOut)
{
    return TwoSpinIn[0] * SPIN3 + TwoSpinIn[1] * SPIN2 +
           TwoSpinOut[0] * SPIN + TwoSpinOut[1];
}

vector<int> Basic::GetSpinIndexVector_FourSpinsFileter(SpinFilter filter)
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

int Basic::TauToBin(real tau)
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

int Basic::TauToBin(real t_in, real t_out)
{
    return TauToBin(t_out - t_in);
}

real Basic::BinToTau(int Bin)
{
    //TODO: mapping between tau and bin
    return Bin * _dBeta + _Beta / 2;
}

void Basic::Reset(real beta)
{
    _Beta = beta;
    _dBeta = beta / MAX_TAU_BIN;
    _dBetaInverse = 1.0 / _dBeta;
    //TODO: please implement how to reset the weight here
}
