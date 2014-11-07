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
            _Weight(i * _Shape[TAU] + j) *= PhaseFactor[j];
}

void WeightNoMeasure::_FFT(fft::Dir direction)
{
    for (int sp = 0; sp < _Shape[SP]; sp++)
        for (int sub = 0; sub < _Shape[SUB]; sub++)
            fft::fftnD((Complex *)_Weight[sp][sub](), (int *)_SpaceTimeShape, _Lat.Dimension + 1, direction);
}

void Sigma::FFT(fft::Dir direction)
{
    if (direction == fft::FORTH)
        _ChangeSymmetry(direction);

    _FFT(direction);

    if (direction == fft::BACK)
        _ChangeSymmetry(direction);
}

void Polar::FFT(fft::Dir direction)
{
    _FFT(direction);
}

void G::FFT(fft::Dir direction)
{
    if (direction == fft::FORTH)
        _ChangeSymmetry(direction);

    _FFT(direction);

    if (direction == fft::BACK)
        _ChangeSymmetry(direction);
}

void W::FFT(fft::Dir direction)
{
    _FFT(direction);
}