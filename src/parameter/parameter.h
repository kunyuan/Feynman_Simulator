//
//  state.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/2/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__state__
#define __Feynman_Simulator__state__

#include "../lattice/lattice.h"
#include "parameter_map.h"
#include "../utility/rng.h"
#include <string>

class Parameter {
  public:
    int PID;
    real Jcp;
    real InitialBeta;
    real DeltaBeta;
    real FinalBeta;
    int Order;

    //derived
    real Beta;
    real T;
    Lattice Lat;

    bool BuildNew(const std::string &InputFile);
    bool Load(const std::string &InputFile);
    void Save(const std::string &InputFile, std::string Mode = "a");

  protected:
    ParameterMap _para;
    Vec<int> L;
};

class ParameterMC : public Parameter {
  public:
    int Version;
    long long Counter;
    int Toss;
    int Sample;
    int Sweep;
    int Seed;
    RandomFactory RNG;
    real WormSpaceReweight;
    real OrderReWeight[MAX_ORDER];

    bool BuildNew(const std::string &InputFile);
    bool Load(const std::string &InputFile);
    void Save(const std::string &InputFile, std::string Mode = "a");
};

#endif /* defined(__Feynman_Simulator__state__) */
