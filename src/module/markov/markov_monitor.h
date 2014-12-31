//
//  measure.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/19/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__measure__
#define __Feynman_Simulator__measure__

#include "module/weight/weight.h"
#include "estimator/estimator.h"
namespace diag {
class Diagram;
}
namespace para {
class ParaMC;
}

namespace mc {
class MarkovMonitor {
public:
    MarkovMonitor();

    para::ParaMC* Para;
    diag::Diagram* Diag;
    weight::Weight* Weight;

    EstimatorBundle<Complex> cEstimator;
    EstimatorBundle<real> rEstimator;
    EstimatorBundle<real> DetailBalanceEstimator;
    Estimator<real> ZeroOrderWeight;

    bool BuildNew(para::ParaMC&, diag::Diagram&, weight::Weight&);
    bool Load(const std::string& InputFile, para::ParaMC&, diag::Diagram&, weight::Weight&);
    void Save(const std::string& InputFile, const std::string& Mode = "a");
    bool FromDict(const Dictionary&, para::ParaMC&, diag::Diagram&, weight::Weight&);
    Dictionary ToDict();
    void ReWeight();

    void Annealing();
    void SqueezeStatistics();
    void ReWeightEachOrder();
    void Measure();
    void AddStatistics();
};
}

#endif /* defined(__Feynman_Simulator__measure__) */
