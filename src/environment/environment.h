//
//  environment.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/17/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__environment__
#define __Feynman_Simulator__environment__

#include "../diagram/diagram.h"
#include "../utility/rng.h"
#include "../observable/weight.h"
#include "../utility/convention.h"
#include "../job/job.h"
#include "../lattice/lattice.h"

class Environment {
  public:
    Environment(Jobs &);
    Lattice Lat;
    real Beta;
    int Order;
    std::string StateFile;

    bool LoadState();
    void SaveState();
};

class EnvMonteCarlo : public Environment {
  public:
    EnvMonteCarlo(Jobs &);
    Diagram Diag;
    long long Counter;
    real OrderWeight[MAX_ORDER];
    EstimatorBundle<Complex> cEstimator;
    EstimatorBundle<real> rEstimator;
    Weight::Sigma Sigma;
    Weight::Polar Polar;
    Weight::W W;
    Weight::G G;
    Weight::Worm Worm;

    void Build();
    bool LoadState();
    void SaveState();
};

class EnvDyson : public Environment {
  public:
    EnvDyson(Jobs &);
};

#endif /* defined(__Feynman_Simulator__environment__) */
