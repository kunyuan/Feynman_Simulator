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

#ifndef GET
#define GET(para, thing)                             \
    {                                                \
        stringstream ss(para.front());               \
        (ss) >> thing;                               \
        if (!(ss).good())                            \
            ABORT("Fail to read " << #thing << "!"); \
        para.pop_front();                            \
    }
#endif

class Environment {
  protected:
    std::list<std::string> _para;
    Vec<int> _L;

  public:
    Environment();
    ~Environment();
    JobType Type;
    int PID;
    real Jcp;
    real InitialBeta;
    real DeltaBeta;
    real Beta;
    int Order;
    bool DoesLoad;
    std::string StateFile;
    Lattice *Lat;

    bool BuildFromFile(std::string InputFile);
    void BuildTest();
    bool LoadState();
    void SaveState();
};

class EnvMonteCarlo : public Environment {
  public:
    EnvMonteCarlo();
    ~EnvMonteCarlo();
    long long Counter;
    int Toss;
    int Sample;
    int Sweep;
    int Seed;
    real WormSpaceReweight;

    Diagram Diag;
    real *OrderWeight;

    EstimatorBundle<Complex> cEstimator;
    EstimatorBundle<real> rEstimator;

    Weight::Worm *WormWeight;
    Weight::Sigma *Sigma;
    Weight::Polar *Polar;
    Weight::W *W;
    Weight::G *G;

    bool BuildFromFile(std::string InputFile);
    bool LoadState();
    void SaveState();
};

class EnvDyson : public Environment {
  public:
    EnvDyson();
    bool BuildFromFile(std::string InputFile);
};

#endif /* defined(__Feynman_Simulator__environment__) */
