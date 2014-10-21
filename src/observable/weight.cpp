//
//  observable.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/13/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "weight.h"

using namespace std;
using namespace Array;
using namespace Weight;

Base::Base(const Lattice &lat, real beta, int order)
{
    _Lat = lat;
    _Beta = beta;
    _Order = order;
}

int Base::SpinIndex(spin SpinIn, spin SpinOut)
{
    return SpinIn * SPIN + SpinOut;
}

int Base::SpinIndex(spin *TwoSpinIn, spin *TwoSpinOut)
{
    return TwoSpinIn[0] * SPIN3 + TwoSpinIn[1] * SPIN2 +
           TwoSpinOut[0] * SPIN + TwoSpinOut[1];
}

int Base::TauToBin(real tau)
{
    //TODO: mapping between tau and bin
    return int(tau * MAX_BIN / _Beta);
}

real Base::BinToTau(int Bin)
{
    //TODO: mapping between tau and bin
    return (real(Bin / MAX_BIN) + 0.5) * _Beta;
}

Sigma::Sigma(const Lattice &lat, real beta, int order)
    : Base(lat, beta, order)
{
    _Shape[ORDER] = _Order;
    _Shape[SP] = SPIN2;
    _Shape[SUB] = lat.SublatVol * lat.SublatVol;
    _Shape[VOL] = lat.Vol;
    _Shape[TAU] = MAX_BIN;

    _Weight = new array4<Complex>(_Shape[SP], _Shape[SUB], _Shape[VOL], _Shape[TAU]);
    _WeightAccu = new array5<Complex>(_Shape[ORDER], _Shape[SP], _Shape[SUB], _Shape[VOL], _Shape[TAU]);
    _WeightSquareAccu = new array2<Complex>(_Shape[ORDER], _Shape[TAU]);

    _Name = "Sigma";
    _Norm = 1.0;
}

Sigma::~Sigma()
{
    delete _Weight;
    delete _WeightAccu;
    delete _WeightSquareAccu;
}
/**
*  \brief Check statistics from StartFromOrder up to _Order
*
*  @param StartFromOrder defines the start order to check
*  @param ErrorThreshold threshold*100% defines how big error is acceptiable
*
*  @return return the maxium order with acceptable errors
*/
int Sigma::OrderAcceptable(int StartFromOrder, real ErrorThreshold)
{
    Complex mean(0.0, 0.0);
    Complex std(0.0, 0.0);
    int order = StartFromOrder;
    for (bool flag = true; flag && order <= _Order; order++) {
        array1<Complex> pSquareAccu = (*_WeightSquareAccu)[order - 1];
        array1<Complex> pAccu = (*_WeightAccu)[order - 1][0][0][0];
        for (int t = 0; t < _Shape[TAU]; t++) {
            mean = pAccu[t] / _Norm;
            std = pSquareAccu[t] / _Norm;
            std = std - mean * mean;
            if (mean.Re > ErrorThreshold * sqrt(std.Re) ||
                mean.Im > ErrorThreshold * sqrt(std.Im)) {
                flag = false;
                break;
            }
        }
    }
    return order - 1;
}

/**
*  Update the Sigma weight up to the UpToOrder
*
*  @param UpToOrder the upper limit of orders to accept
*/

void Sigma::UpdateWeight(int UpToOrder)
{
    int size = _Weight->Size();
    int order = 1;
    for (int i = 0; i < size; i++)
        //assign order=1 directly to initialize _Weight
        (*_Weight)(i) = (*_WeightAccu)[order - 1](i) / _Norm;

    for (order = 2; order <= UpToOrder; order++) {
        //add order>1 on _Weight
        for (int i = 0; i < size; i++)
            (*_Weight)(i) += (*_WeightAccu)[order - 1](i) / _Norm;
    }
}

Complex Sigma::Weight(const Distance &d, real dtau, spin SpinIn, spin SpinOut)
{
    return (*_Weight)[SpinIndex(SpinIn, SpinOut)][d.SublatIndex][d.CoordiIndex][TauToBin(dtau)];
}

//Estimate<Complex> Sigma::WeightWithError(const Distance &d, real dtau, spin SpinIn, spin SpinOut)
//{
//}

void Sigma::Measure(const Distance &d, real dtau, spin SpinIn, spin SpinOut, int order, const Complex &weight)
{
    if (DEBUGMODE && order <= 0)
        LOG_ERROR("Too small order=" << order << endl);
    int spin_index = SpinIndex(SpinIn, SpinOut);
    int tau_bin = TauToBin(dtau);
    (*_WeightAccu)[order - 1][spin_index][d.SublatIndex][d.CoordiIndex][tau_bin] += weight;
    if (spin_index == 0 && d.SublatIndex == 0 && d.CoordiIndex == 0)
        (*_WeightSquareAccu)[order - 1][tau_bin] += weight * weight;
    _Norm += 1.0;
}

void Sigma::ClearStatistics()
{
    _Norm = 1.0;
    int size = _WeightAccu->Size();
    for (int i = 0; i < size; i++)
        (*_WeightAccu)(i) = 0.0;
    size = _WeightSquareAccu->Size();
    for (int i = 0; i < size; i++)
        (*_WeightSquareAccu)(i) = 0.0;
}

void Sigma::SqueezeStatistics(real factor)
{
    _Norm /= factor;
    int size = _WeightAccu->Size();
    for (int i = 0; i < size; i++)
        (*_WeightAccu)(i) /= factor;
    size = _WeightSquareAccu->Size();
    for (int i = 0; i < size; i++)
        (*_WeightSquareAccu)(i) /= factor;
}

/************************   Polarization   *********************************/
//
Polar::Polar(const Lattice &lat, real beta, int order)
    : Base(lat, beta, order)
{
    _Shape[ORDER] = _Order;
    _Shape[SP] = SPIN4;
    _Shape[SUB] = lat.SublatVol * lat.SublatVol;
    _Shape[VOL] = lat.Vol;
    _Shape[TAU] = MAX_BIN;

    _Weight = new array4<Complex>(_Shape[SP], _Shape[SUB], _Shape[VOL], _Shape[TAU]);
    _WeightAccu = new array5<Complex>(_Shape[ORDER], _Shape[SP], _Shape[SUB], _Shape[VOL], _Shape[TAU]);
    _WeightSquareAccu = new array2<Complex>(_Shape[ORDER], _Shape[TAU]);

    _Name = "Polar";
    _Norm = 1.0;
}

Polar::~Polar()
{
    delete _Weight;
    delete _WeightAccu;
    delete _WeightSquareAccu;
}

//G::G(const Lattice &lat, real beta)
//    : Base(lat, beta)
//{
//    _Weight = new Array4<Complex>(SPIN2, lat.SublatVol * lat.SublatVol, lat.Vol, MAX_BIN);
//}
//
//G::~G()
//{
//    delete _Weight;
//}
//
//Complex G::Weight(const Distance &dR, real dtau, spin SpinIn, spin SpinOut)
//{
//    return (*_Weight)[SpinIndex(SpinIn, SpinOut)][dR.SublatIndex][dR.CoordiIndex][TauToBin(dtau)];
//}
//
//W::W(const Lattice &lat, real beta)
//    : Base(lat, beta)
//{
//    _Weight = new Array4<Complex>(SPIN4, lat.SublatVol * lat.SublatVol, lat.Vol, MAX_BIN);
//}
//
//Weight::W::~W()
//{
//    delete _Weight;
//}
