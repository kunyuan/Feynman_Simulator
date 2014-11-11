//
//  weightEstimator.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/10/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "weight_estimator.h"
#include "utility/abort.h"
#include "utility/scopeguard.h"

using namespace std;
using namespace Array;
using namespace weight;

/**********************   Weight Needs measuring  **************************/

WeightNeedMeasure::WeightNeedMeasure(const Lattice &lat, real beta, int order,
                                     bool IsTauSymmetric, int SpinVol,
                                     string name, real Norm)
    : WeightNoMeasure(lat, beta, order, IsTauSymmetric, SpinVol, name)
{
    _Norm = Norm;
    //use _Shape[ORDER] to _Shape[TAU] to construct array5
    _WeightAccu.Allocate(Shape());
    for (int i = 1; i <= order; i++)
        _WeightErrorEstimator.AddEstimator(name + "_AvgofOrder" + ToString(i));
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

Complex WeightNeedMeasure::RelativeError(int order)
{
    return _WeightErrorEstimator[order - 1].Estimate().RelativeError();
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
    int order = StartFromOrder;
    for (; order <= _Order; order++) {
        Complex rerr = _WeightErrorEstimator[order - 1].Estimate().RelativeError();
        if (rerr.Re > ErrorThreshold || rerr.Im > ErrorThreshold)
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
    _WeightErrorEstimator.AddStatistics();
}

void WeightNeedMeasure::ClearStatistics()
{
    _NormAccu = 1.0;
    _WeightAccu = 0.0;
    _WeightErrorEstimator.ClearStatistics();
}
//TODO: you may have to replace int with size_t here

void WeightNeedMeasure::SqueezeStatistics(real factor)
{
    if (DEBUGMODE && factor <= 0.0)
        ABORT("factor=" << factor << "<=0!");
    _NormAccu /= factor;
    _WeightAccu *= 1.0 / factor;
    _WeightErrorEstimator.SqueezeStatistics(factor);
}

/**********************   Weight IO ****************************************/
void WeightNeedMeasure::Save(const std::string &FileName, std::string Mode)
{
    unsigned int shape[1] = {1};
    cnpy::npz_save(FileName, _Name + ".Norm", &_Norm, shape, 1, Mode);
    cnpy::npz_save(FileName, _Name + ".NormAccu", &_NormAccu, shape, 1, "a");
    cnpy::npz_save(FileName, _Name + ".WeightAccu", _WeightAccu(), Shape(), 5, "a");
    WeightNoMeasure::Save(FileName);
    _WeightErrorEstimator.SaveStatistics(FileName, "a");
}

bool WeightNeedMeasure::Load(const std::string &FileName)
{
    _WeightErrorEstimator.LoadStatistics(FileName);
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