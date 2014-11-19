//
//  weightEstimator.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/10/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__weightEstimator__
#define __Feynman_Simulator__weightEstimator__

#include "weight_base.h"
#include "estimator/estimator.h"

namespace weight {

class WeightNeedMeasure : public WeightNoMeasure {
  public:
    WeightNeedMeasure(const Lattice &, real Beta, int order, bool IsTauSymmetric,
                      int Spine, std::string, real Norm);

    //Shape including ORDER
    unsigned int *ExtendedShape();

    //The normalization facto is not considered in _WeightErrorEstimator,
    //thus only relative error makes sense
    Complex RelativeError(int order);

    int OrderAcceptable(int StartFromOrder, real ErrorThreshold);
    //update final weight density to WeightNoMeasure._Weight
    void UpdateWeight(int UpToOrder);

    //The internal _Beta will be changed, so do _WeightAccu, _DeltaWeightAccu and _NormAccu
    //all changed will be done to make sure GetWeightArray returns the reweighted weight function
    //(as for now, reweighted weight function is set to be the unreweighted weight function)
    void ReWeight(real Beta);

    void MeasureNorm();

    //add statistics to the history of _WeightErrorEstimator, so that
    //weight error can be estimated. !!!It cosumes memory!!!
    void AddStatistics();
    void ClearStatistics();
    void SqueezeStatistics(real factor);
    //    std::string PrettyString();
    void Save(const std::string &FileName, std::string Mode = "a");
    bool Load(const std::string &);

  protected:
    unsigned int _MeaShape[5];
    int _Order;
    real _Norm;     //The normalization factor
    real _NormAccu; //The normalization accumulation
    //final weight of each bin = _WeightAccu/_NormAccu*_Norm
    //final weight function = (final weight of each bin)*MAX_BIN/Beta
    Array::array5<Complex> _WeightAccu; //dim=0 is order, than follows WeightNoMeasure::Shape()
    EstimatorBundle<Complex> _WeightErrorEstimator;
};
}
#endif /* defined(__Feynman_Simulator__weightEstimator__) */
