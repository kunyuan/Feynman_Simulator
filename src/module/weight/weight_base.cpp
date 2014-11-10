//
//  observable.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/13/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "weight.h"
#include "utility/abort.h"
#include "utility/scopeguard.h"

using namespace std;
using namespace Array;
using namespace weight;

WeightNoMeasure::WeightNoMeasure(const Lattice &lat, real beta,
                                 int order, int SpinVol, string name)
    : _Lat(lat),
      _Beta(beta),
      _Order(order),
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

    BareWeight.Allocate(Shape());
    BareWeight = 0.0;
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

int WeightNoMeasure::SpinIndex(spin *TwoSpinIn, spin *TwoSpinOut)
{
    return TwoSpinIn[0] * SPIN3 + TwoSpinIn[1] * SPIN2 +
           TwoSpinOut[0] * SPIN + TwoSpinOut[1];
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

int WeightNoMeasure::TauToBin(real t_in, real t_out)
{
    return TauToBin(t_out - t_in);
}

real WeightNoMeasure::BinToTau(int Bin)
{
    //TODO: mapping between tau and bin
    return Bin * _dBeta + _Beta / 2;
}

void WeightNoMeasure::SetTest()
{
    for (unsigned int i = 0; i < SmoothWeight.Size(); i++) {
        SmoothWeight(i) = Complex(1.0, 0.0);
    }
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

/**********************   Weight Needs measuring  **************************/

WeightNeedMeasure::WeightNeedMeasure(const Lattice &lat, real beta,
                                     int order, int SpinVol, string name, real Norm)
    : WeightNoMeasure(lat, beta, order, SpinVol, name)
{
    _Norm = Norm;
    //use _Shape[ORDER] to _Shape[TAU] to construct array5
    _WeightAccu.Allocate(Shape());
    for (int i = 1; i <= order; i++)
        _Average.AddEstimator(name + "_AvgofOrder" + ToString(i));
    ClearStatistics();
}

unsigned int *WeightNeedMeasure::Shape()
{
    return _Shape;
}

void WeightNeedMeasure::ReWeight(real Beta)
{
    //make sure
    //real NormFactor = 1.0 / _NormAccu * _Norm * MAX_BIN / _Beta;
    //has the same value before Beta is changed
    //so that GetWeightArray will give a same weight function
    _NormAccu *= _Beta / Beta;
    _Beta = Beta;
    _dBeta = Beta / MAX_BIN;
    _dBetaInverse = 1.0 / _dBeta;
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
    //WeightEstimator and the corresponding WeightArray
    //share the same memory structure from _Shape[SP]  to _Shape[TAU]
    int order = 1;
    real NormFactor = 1.0 / _NormAccu * _Norm * MAX_BIN / _Beta;

    SmoothWeight = _WeightAccu[order - 1];
    //add order>1 on _Weight
    for (order = 2; order <= UpToOrder; order++)
        SmoothWeight += _WeightAccu[order - 1];
    SmoothWeight *= NormFactor;
}

void WeightNeedMeasure::MeasureNorm()
{
    _NormAccu += 1.0;
}

void WeightNeedMeasure::AddStatistics()
{
    _Average.AddStatistics();
}

void WeightNeedMeasure::ClearStatistics()
{
    _NormAccu = 1.0;
    _WeightAccu = 0.0;
    _Average.ClearStatistics();
}
//TODO: you may have to replace int with size_t here

void WeightNeedMeasure::SqueezeStatistics(real factor)
{
    if (DEBUGMODE && factor <= 0.0)
        ABORT("factor=" << factor << "<=0!");
    _NormAccu /= factor;
    _WeightAccu *= 1.0 / factor;
    _Average.SqueezeStatistics(factor);
}

/**********************   Weight IO ****************************************/
void WeightNeedMeasure::Save(const std::string &FileName, std::string Mode)
{
    unsigned int shape[1] = {1};
    cnpy::npz_save(FileName, _Name + ".Norm", &_Norm, shape, 1, Mode);
    cnpy::npz_save(FileName, _Name + ".NormAccu", &_NormAccu, shape, 1, "a");
    cnpy::npz_save(FileName, _Name + ".WeightAccu", _WeightAccu(), Shape(), 5, "a");
    WeightNoMeasure::Save(FileName);
    _Average.SaveStatistics(FileName, "a");
}

bool WeightNeedMeasure::Load(const std::string &FileName)
{
    _Average.LoadStatistics(FileName);
    WeightNoMeasure::Load(FileName);

    cnpy::npz_t NpzMap = cnpy::npz_load(FileName);
    ON_SCOPE_EXIT([&] {NpzMap.destruct(); });

    cnpy::NpyArray weight_accu = NpzMap[_Name + ".WeightAccu"];
    if (weight_accu.data == nullptr)
        ABORT("Can't find estimator " << _Name << ".WeightAccu in .npz data file!");
    //using assign here will make a copy of the data in Complex *start
    _WeightAccu = reinterpret_cast<Complex *>(weight_accu.data);

    //read normalization factor
    cnpy::npz_load_number(NpzMap, _Name + ".NormAccu", _NormAccu);
    cnpy::npz_load_number(NpzMap, _Name + ".Norm", _Norm);

    return true;
}
