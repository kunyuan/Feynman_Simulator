//
//  weight_estimator.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/25/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "utility/abort.h"
#include "utility/scopeguard.h"
#include "utility/dictionary.h"
#include "weight_estimator.h"

using namespace std;
using namespace Array;
using namespace weight;

/**********************   Weight Needs measuring  **************************/

WeightEstimator::WeightEstimator(const Lattice& lat, real beta, int order, string name, real Norm, const uint* Shape)
{
    _Vol = lat.Vol;
    _SublatVol = lat.SublatVol;
    _MaxTauBin = Shape[TAU];
    _Beta = beta;
    _dBeta = beta / _MaxTauBin;
    _dBetaInverse = 1.0 / _dBeta;
    _Order = order;
    _Name = name;
    _MeaShape[0] = _Order;
    std::copy(Shape, Shape + 4, &_MeaShape[1]);
    //use _MeaShape[0] to _MeaShape[TAU] to construct array5
    _WeightAccu.Allocate(_MeaShape);
    _Norm = Norm * (_MaxTauBin / _Beta) / _Beta / _Vol / _SublatVol;
    ClearStatistics();
}

void WeightEstimator::Anneal(real Beta)
{
    //make sure
    //real NormFactor = 1.0 / _NormAccu * _Norm;
    //has the same value before Beta is changed
    //so that GetWeightArray will give a same weight function
    _NormAccu *= pow((Beta / _Beta), 2.0);
}

void WeightEstimator::MeasureNorm()
{
    _NormAccu += 1.0;
}

void WeightEstimator::Measure(uint* Index, int Order, Complex weight)
{
    if (DEBUGMODE && Order < 1)
        LOG_ERROR("Too small order=" << Order);
    _WeightAccu[Order - 1][Index[SP]][Index[SUB]]
               [Index[VOL]][Index[TAU]] += weight;
}

void WeightEstimator::ClearStatistics()
{
    _NormAccu = 0.0;
    _WeightAccu = 0.0;
}
//TODO: you may have to replace int with size_t here

void WeightEstimator::SqueezeStatistics(real factor)
{
    if (DEBUGMODE && factor <= 0.0)
        ABORT("factor=" << factor << "<=0!");
    _NormAccu /= factor;
    _WeightAccu *= 1.0 / factor;
}

/**********************   Weight IO ****************************************/

bool WeightEstimator::FromDict(const Dictionary& dict)
{
    _Norm = dict.Get<real>("Norm");
    _NormAccu = dict.Get<real>("NormAccu");
    auto arr = dict.Get<Python::ArrayObject>("WeightAccu");
    ASSERT_ALLWAYS(Equal(arr.Shape().data(), _MeaShape, 5), "Shape should match!");
    _WeightAccu = arr.Data<Complex>();
    return true;
}

Dictionary WeightEstimator::ToDict()
{
    Dictionary dict;
    dict["Norm"] = _Norm;
    dict["NormAccu"] = _NormAccu;
    dict["WeightAccu"] = Python::ArrayObject(_WeightAccu(), _MeaShape, 5);
    return dict;
}
