//
//  weight_estimator.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/25/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__weight_estimator__
#define __Feynman_Simulator__weight_estimator__

#include "estimator/estimator.h"
#include "index_map.h"
#include "weight_array.h"

class Dictionary;
namespace weight {

class IndexMap;
class WeightEstimator {
public:
    WeightEstimator();
    void Allocate(const IndexMap& map, int order, real Norm);

    //The internal _Beta will be changed, so do _WeightAccu, _DeltaWeightAccu and _NormAccu
    //all changed will be done to make sure GetWeightArray returns the reweighted weight function
    //(as for now, reweighted weight function is set to be the unreweighted weight function)
    void Anneal(real Beta);

    void MeasureNorm(real weight);
    void Measure(uint WeightIndex, int Order, Complex Weight);

    void ClearStatistics();
    void SqueezeStatistics(real factor);
    //    std::string PrettyString();
    bool FromDict(const Dictionary&);
    Dictionary ToDict();

protected:
    real _Beta;
    real _Norm; //The normalization factor
    real _NormAccu; //The normalization accumulation
    //final weight function =_WeightAccu/_NormAccu*_Norm
    //final weight of each bin = (final weight of each bin)/MAX_BIN*Beta
    WeightArray<SMOOTH_T_SIZE + 1> _WeightAccu; //dim=0 is order
    uint _WeightSize;
};
}
#endif /* defined(__Feynman_Simulator__weight_estimator__) */
