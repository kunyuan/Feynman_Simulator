//
//  environment.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/17/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__environment__
#define __Feynman_Simulator__environment__

#include "job.h"
#include "lattice.h"
#include "diagram.h"
#include "rng.h"
#include "observable.h"

class Environment
{
public:
    Environment(Jobs &);
    real Beta;
    int Order;
    Lattice Lat;
    Diagram Diag;
    RandomFactory *RNG;
    
    EstimatorBundle<Complex> cEstimator;
    EstimatorBundle<real> rEstimator;
    
    void Build();
    bool ReadState();
    void WriteState();
};

class EnvMoneCarlo: public Environment
{
public:
    EnvMoneCarlo(Jobs &);
    real *OrderWeight;
};

class EnvDyson: public Environment
{
public:
    EnvDyson(Jobs &);
};

#endif /* defined(__Feynman_Simulator__environment__) */
