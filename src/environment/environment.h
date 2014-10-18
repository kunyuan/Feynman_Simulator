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
#include "diagram_object.h"
#include "convention.h"

class Environment
{
public:
    Environment(Jobs &);
    real Beta;
    int Order;
    Lattice Lat;
    
    bool ReadState();
    void SaveState();
};

class EnvMoneCarlo: public Environment
{
public:
    EnvMoneCarlo(Jobs &);
    RandomFactory RNG;
    Diagram Diag;
    int Counter;
    real OrderWeight[MAX_ORDER];
    EstimatorBundle<Complex> cEstimator;
    EstimatorBundle<real> rEstimator;
    
    void Annealing();
    void SqueezeStatistics();
    void ReWeightEachOrder();
    void Measure();
    void AddStatistics();
    
    void Build();
    bool ReadState();
    void SaveState();
};

class EnvDyson: public Environment
{
public:
    EnvDyson(Jobs &);
};

#endif /* defined(__Feynman_Simulator__environment__) */
