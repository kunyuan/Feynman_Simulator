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

class Dictionary;
namespace weight {

class WeightEstimator {
public:
    WeightEstimator(const Lattice& lat, real beta, int order, string name, real Norm, const uint* Shape);
    //Shape including ORDER
    uint* GetExtendedShape();

    //The internal _Beta will be changed, so do _WeightAccu, _DeltaWeightAccu and _NormAccu
    //all changed will be done to make sure GetWeightArray returns the reweighted weight function
    //(as for now, reweighted weight function is set to be the unreweighted weight function)
    void Anneal(real Beta);

    void MeasureNorm();
    void Measure(uint* Index, int Order, Complex Weight);

    //add statistics to the history of _WeightErrorEstimator, so that
    //weight error can be estimated. !!!It cosumes memory!!!
    void AddStatistics();
    void ClearStatistics();
    void SqueezeStatistics(real factor);
    //    std::string PrettyString();
    bool FromDict(const Dictionary&);
    Dictionary ToDict();

protected:
    unsigned int _MeaShape[5];
    uint _MaxTauBin;
    int _Vol;
    int _SublatVol;
    real _Beta;
    real _dBeta; //_Beta/MAX_TAU
    real _dBetaInverse;
    int _Order;
    real _Norm; //The normalization factor
    real _NormAccu; //The normalization accumulation
    //final weight function =_WeightAccu/_NormAccu*_Norm
    //final weight of each bin = (final weight of each bin)/MAX_BIN*Beta
    Array::array5<Complex> _WeightAccu; //dim=0 is order, than follows WeightNoMeasure::Shape()
    std::string _Name;
};
}
#endif /* defined(__Feynman_Simulator__weight_estimator__) */
