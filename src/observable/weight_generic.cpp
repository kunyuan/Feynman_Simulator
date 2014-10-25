//
//  observable.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/13/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "weight.h"
#include "../utility/utility.h"
#include "../utility/abort.h"

using namespace std;
using namespace Array;
using namespace Weight;

WeightNoMeasure::WeightNoMeasure(const Lattice &lat, real beta, int order, int SpinVol, string name)
{
    _Lat = lat;
    _Beta = beta;
    _Order = order;
    _dBeta = beta / MAX_BIN;
    _dBetaInverse = MAX_BIN / beta;
    _Name = name;

    _Shape[ORDER] = order;
    _Shape[SP] = SpinVol;
    _Shape[SUB] = lat.SublatVol * lat.SublatVol;
    _Shape[VOL] = lat.Vol;
    _Shape[TAU] = MAX_BIN;

    _Weight.Allocate((unsigned int *)(_Shape + SP));
    //use _Shape[SP] to _Shape[TAU] to construct array4
}

int WeightNoMeasure::SpinIndex(spin SpinIn, spin SpinOut)
{
    return SpinIn * SPIN + SpinOut;
}

int WeightNoMeasure::SpinIndex(spin *TwoSpinIn, spin *TwoSpinOut)
{
    return TwoSpinIn[0] * SPIN3 + TwoSpinIn[1] * SPIN2 +
           TwoSpinOut[0] * SPIN + TwoSpinOut[1];
}

int WeightNoMeasure::TauToBin(real tau)
{
    //TODO: mapping between tau and bin
    if (tau < 0) {
        if (DEBUGMODE && tau < -_Beta)
            LOG_ERROR("Beta=" << tau << " is too small!");
        return int(tau * _dBetaInverse) + MAX_BIN;
    }
    else {
        if (DEBUGMODE && tau >= _Beta)
            LOG_ERROR("Beta=" << tau << " is too large!");
        return int(tau * _dBetaInverse);
    }
}

real WeightNoMeasure::BinToTau(int Bin)
{
    //TODO: mapping between tau and bin
    return Bin * _dBeta + _Beta / 2;
}

void WeightNoMeasure::InitializeState()
{
    for (unsigned int i = 0; i < _Weight.Size(); i++) {
        _Weight(i) = Complex(2.0, 0.0);
    }
}

void WeightNoMeasure::SaveState(const std::string &FileName, std::string Mode)
{
    cnpy::npz_save(cnpy::npz_name(FileName), _Name, _Weight.Data(), _Shape + SP, 4, Mode);
}

bool WeightNoMeasure::LoadState(const std::string &FileName)
{
    cnpy::NpyArray weight = cnpy::npz_load(cnpy::npz_name(FileName), _Name);
    if (weight.data == nullptr)
        ABORT("Can't find estimator " << _Name << " in .npz data file!" << endl);
    _Weight = reinterpret_cast<Complex *>(weight.data);
    return true;
}

/**********************   Weight Needs measuring  **************************/

WeightNeedMeasure::WeightNeedMeasure(const Lattice &lat, real beta, int order, int SpinVol, string name)
    : WeightNoMeasure(lat, beta, order, SpinVol, name)
{
    _WeightAccu.Allocate((unsigned int *)_Shape);
    //use _Shape[ORDER] to _Shape[TAU] to construct array5
    _Norm = _dBeta;
    for (int i = 1; i <= order; i++)
        _Average.AddEstimator(name + "_AvgofOrder" + ToString(i));
}

Estimate<Complex> WeightNeedMeasure::WeightWithError(int order)
{
    return _Average[order - 1].Estimate();
}

/**
*  \brief Check statistics from StartFromOrder up to _Order
*
*  @param StartFromOrder defines the start order to check
*  @param ErrorThreshold threshold*100% defines how big error is acceptiable
*
*  @return return the maxium order with acceptable errors
*/
int WeightNeedMeasure::OrderAcceptable(int StartFromOrder, real ErrorThreshold)
{
    bool flag = true;
    int order = StartFromOrder;
    Estimate<Complex> temp;
    for (; order <= _Order; order++) {
        temp = _Average[order - 1].Estimate();
        flag = (temp.Mean.Re > temp.Error.Re * ErrorThreshold) &&
               (temp.Mean.Im > temp.Error.Im * ErrorThreshold);
        if (!flag)
            break;
    }
    return order - 1;
}

/**
*  Update the Sigma weight up to the UpToOrder
*
*  @param UpToOrder the upper limit of orders to accept
*/

void WeightNeedMeasure::UpdateWeight(int UpToOrder)
{
    int size = _Weight.Size();
    int order = 1;
    for (int i = 0; i < size; i++)
        //assign order=1 directly to initialize _Weight
        _Weight(i) = _WeightAccu[order - 1](i) / _Norm;

    for (order = 2; order <= UpToOrder; order++) {
        //add order>1 on _Weight
        for (int i = 0; i < size; i++)
            _Weight(i) += _WeightAccu[order - 1](i) / _Norm;
    }
}

void WeightNeedMeasure::AddStatistics()
{
    _Average.AddStatistics();
}

void WeightNeedMeasure::ClearStatistics()
{
    _Norm = _dBeta;
    int size = _WeightAccu.Size();
    for (int i = 0; i < size; i++)
        _WeightAccu(i) = 0.0;
    _Average.ClearStatistics();
}

void WeightNeedMeasure::SqueezeStatistics(real factor)
{
    if (DEBUGMODE && factor <= 0.0)
        ABORT("factor=" << factor << "<=0!" << endl);
    _Norm /= factor;
    int size = _WeightAccu.Size();
    for (int i = 0; i < size; i++)
        _WeightAccu(i) /= factor;
    _Average.SqueezeStatistics(factor);
}

/**********************   Weight IO ****************************************/
void WeightNeedMeasure::SaveState(const std::string &FileName, std::string Mode)
{
    unsigned int shape[1] = {1};
    cnpy::npz_save(cnpy::npz_name(FileName), _Name + "_Norm", &_Norm, shape, 1, Mode);
    cnpy::npz_save(cnpy::npz_name(FileName), _Name + "_Accu", _WeightAccu(), _Shape, 5, "a");
    WeightNoMeasure::SaveState(FileName);
    _Average.SaveState(FileName, "a");
}

bool WeightNeedMeasure::LoadState(const std::string &FileName)
{
    _Average.LoadState(FileName);
    WeightNoMeasure::LoadState(FileName);

    cnpy::npz_t NpzMap = cnpy::npz_load(cnpy::npz_name(FileName));
    cnpy::NpyArray sigma_accu = NpzMap[_Name + "_Accu"];
    if (sigma_accu.data == nullptr)
        ABORT("Can't find estimator " << _Name << " _Accu in .npz data file!" << endl);
    _WeightAccu = reinterpret_cast<Complex *>(sigma_accu.data);
    //using assign here will make a copy of the data in Complex *start

    //read normalization factor
    cnpy::NpyArray norm = NpzMap[_Name + "_Norm"];
    if (norm.data == nullptr)
        ABORT("Can't find estimator " << _Name << "_Norm in .npz data file!" << endl);
    _Norm = *reinterpret_cast<real *>(norm.data);

    NpzMap.destruct();
    return true;
}
