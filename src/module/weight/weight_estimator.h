//
//  weight_estimator.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/25/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__weight_estimator__
#define __Feynman_Simulator__weight_estimator__

#include "weight_basic.h"
#include "estimator/estimator.h"

namespace weight {

class WeightEstimator {
public:
    WeightEstimator(real beta, int order, string name, real Norm, const uint* Shape);
    //Shape including ORDER
    uint* GetExtendedShape();

    //The normalization facto is not considered in _WeightErrorEstimator,
    //thus only relative error makes sense
    Complex RelativeError(int order);

    int OrderAcceptable(int StartFromOrder, real ErrorThreshold);
    //update final weight density to WeightNoMeasure._Weight
    void UpdateWeight(SmoothTMatrix&, int UpToOrder);

    //The internal _Beta will be changed, so do _WeightAccu, _DeltaWeightAccu and _NormAccu
    //all changed will be done to make sure GetWeightArray returns the reweighted weight function
    //(as for now, reweighted weight function is set to be the unreweighted weight function)
    void ReWeight(real Beta);

    void MeasureNorm();
    void Measure(uint* Index, int Order, Complex Weight);

    //add statistics to the history of _WeightErrorEstimator, so that
    //weight error can be estimated. !!!It cosumes memory!!!
    void AddStatistics();
    void ClearStatistics();
    void SqueezeStatistics(real factor);
    //    std::string PrettyString();
    void Save(const std::string& FileName, const std::string& Mode = "a");
    bool Load(const std::string& FileName);

protected:
    unsigned int _MeaShape[5];
    uint _MaxTauBin;
    real _Beta;
    real _dBeta; //_Beta/MAX_TAU
    real _dBetaInverse;
    int _Order;
    real _Norm; //The normalization factor
    real _NormAccu; //The normalization accumulation
    //final weight of each bin = _WeightAccu/_NormAccu*_Norm
    //final weight function = (final weight of each bin)*MAX_BIN/Beta
    Array::array5<Complex> _WeightAccu; //dim=0 is order, than follows WeightNoMeasure::Shape()
    EstimatorBundle<Complex> _WeightErrorEstimator;
    std::string _Name;
};
}
#endif /* defined(__Feynman_Simulator__weight_estimator__) */
