//
//  measure.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/19/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__measure__
#define __Feynman_Simulator__measure__

#include "../observable/weight.h"
#include "../parameter/parameter.h"
#include "../diagram/diagram.h"

class MarkovMonitor {
  public:
    MarkovMonitor();

    Parameter *Para;
    Diagram *Diag;
    weight::Weight *Weight;

    EstimatorBundle<Complex> cEstimator;
    EstimatorBundle<real> rEstimator;
    EstimatorBundle<real> DetailBalanceEstimator;
    Estimator<real> ZeroOrderWeight;

    bool BuildNew(ParameterMC &, Diagram &, weight::Weight &);
    bool Load(const std::string &InputFile, ParameterMC &, Diagram &, weight::Weight &);
    void Save(const std::string &InputFile, const std::string &Mode = "a");

    void Annealing();
    void SqueezeStatistics();
    void ReWeightEachOrder();
    void Measure();
    void AddStatistics();
};
#endif /* defined(__Feynman_Simulator__measure__) */
