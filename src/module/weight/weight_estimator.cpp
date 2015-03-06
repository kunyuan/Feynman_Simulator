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
using namespace weight;

/**********************   Weight Needs measuring  **************************/

WeightEstimator::WeightEstimator()
{
}

void WeightEstimator::Allocate(const IndexMap& map, int order, real Norm)
{
    int Vol = map.Lat.Vol;
    _Beta = map.Beta;
    _Norm = Norm * (map.MaxTauBin / _Beta) / _Beta / Vol;
    uint MeaShape[SMOOTH_T_SIZE + 1];
    MeaShape[0] = order;
    std::copy(map.GetShape(), map.GetShape() + SMOOTH_T_SIZE, &MeaShape[1]);
    _WeightAccu.Allocate(MeaShape, SMOOTH);
    _WeightSize = _WeightAccu.GetSize() / order;
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

void WeightEstimator::MeasureNorm(real weight)
{
    _NormAccu += weight;
}

void WeightEstimator::Measure(uint WeightIndex, int Order, Complex weight)
{
    if (DEBUGMODE && Order < 1)
        LOG_ERROR("Too small order=" << Order);
    uint Index = (Order - 1) * _WeightSize + WeightIndex;
    _WeightAccu[Index] += weight;
}

void WeightEstimator::ClearStatistics()
{
    _NormAccu = 0.0;
    _WeightAccu.Assign(0.0);
}
//TODO: you may have to replace int with size_t here

void WeightEstimator::SqueezeStatistics(real factor)
{
    ASSERT_ALLWAYS(factor > 0, "factor=" << factor << "<=0!");
    _NormAccu /= factor;
    _WeightAccu *= 1.0 / factor;
}

/**********************   Weight IO ****************************************/

bool WeightEstimator::FromDict(const Dictionary& dict)
{
    _Norm = dict.Get<real>("Norm");
    _NormAccu = dict.Get<real>("NormAccu");
    auto arr = dict.Get<Python::ArrayObject>("WeightAccu");
    //assert estimator shape except order dimension
    ASSERT_ALLWAYS(Equal(arr.Shape().data() + 1, _WeightAccu.GetShape() + 1, _WeightAccu.GetDim() - 1), "Shape should match!");
    _WeightAccu.Assign(0.0);
    _WeightAccu.Assign(arr.Data<Complex>(), arr.Size());
    return true;
}

Dictionary WeightEstimator::ToDict()
{
    Dictionary dict;
    dict["Norm"] = _Norm;
    dict["NormAccu"] = _NormAccu;
    dict["WeightAccu"] = Python::ArrayObject(_WeightAccu.Data(), _WeightAccu.GetShape(), _WeightAccu.GetDim());
    return dict;
}
