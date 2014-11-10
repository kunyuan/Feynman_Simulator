//
//  fft.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/5/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "weight.h"

using namespace weight;

int GetSign(fft::Dir direction)
{
    return (direction == fft::FORTH) ? 1 : -1;
}

/**
*  fft::FORTH would change antisymmetric _Weight to be symmetric, and always performed in Tau space
*  fft::BACK would change symmetric _Weight to be antisymmetric, and always performed in Omega space
*/
void WeightNoMeasure::_ChangeSymmetry(fft::Dir direction)
{
    int sign = GetSign(direction);
    Complex PhaseFactor[MAX_BIN];
    for (int i = 0; i < _Shape[TAU]; i++)
        PhaseFactor[i] = exp(Complex(0.0, sign * BinToTau(i) * PI));

    int NumOfTimeSeries = _Shape[SP] * _Shape[SUB] * _Shape[VOL];
    for (int i = 0; i < NumOfTimeSeries; i += _Shape[TAU])
        for (int j = 0; j < _Shape[TAU]; j++)
            SmoothWeight(i * _Shape[TAU] + j) *= PhaseFactor[j];
}

void WeightNoMeasure::_FFT(fft::Dir direction, Mode mode)
{
    //mode==Spatial, all spatial elements of DoIt will be set to be true, so that fft would be performed on it
    //same for mode==Time
    bool DoIt[_Lat.Dimension + 1];
    for (auto i : DoIt)
        DoIt[i] = false;
    if (mode & Spatial)
        for (int i = 0; i < _Lat.Dimension; i++)
            DoIt[i] = true;
    if (mode & Time)
        DoIt[_Lat.Dimension] = true;

    for (int sp = 0; sp < _Shape[SP]; sp++)
        for (int sub = 0; sub < _Shape[SUB]; sub++)
            fft::fftnD((Complex *)SmoothWeight[sp][sub](), (int *)_SpaceTimeShape, _Lat.Dimension + 1, direction, DoIt);
}

void Sigma::FFT(fft::Dir direction, Mode mode)
{
    if (mode & Time && direction == fft::FORTH)
        _ChangeSymmetry(direction);

    _FFT(direction, mode);

    if (mode & Time && direction == fft::BACK)
        _ChangeSymmetry(direction);
}

void Polar::FFT(fft::Dir direction, Mode mode)
{
    _FFT(direction, mode);
}

void G::FFT(fft::Dir direction, Mode mode)
{
    if (mode & Time && direction == fft::FORTH)
        _ChangeSymmetry(direction);

    _FFT(direction, mode);

    if (mode & Time && direction == fft::BACK)
        _ChangeSymmetry(direction);
}

void W::FFT(fft::Dir direction, Mode mode)
{
    _FFT(direction, mode);
}