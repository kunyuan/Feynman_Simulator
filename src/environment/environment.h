//
//  environment.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/17/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__environment__
#define __Feynman_Simulator__environment__

#include "../parameter/parameter.h"
#include "../diagram/diagram.h"
#include "../utility/rng.h"
#include "../observable/weight.h"
#include "../markov/markov_monitor.h"
#include "../markov/markov.h"
#include "../lattice/lattice.h"
#include "../utility/scopeguard.h"

class Environment {
  public:
    int PID;
    enum flag { Bare,
                OldGW };

  protected:
    Environment(int pid);
    std::string _ParaFile();
    std::string _ControlFile();
    std::string _WeightFile();
    std::string _StatisFile();
};

class EnvMonteCarlo : public Environment {
  public:
    EnvMonteCarlo(int pid);
    //can be read from StateFile or InputFile
    ParameterMC Para;
    weight::Weight Weight;
    Diagram Diag;
    Markov Grasshopper;
    MarkovMonitor Scarecrow;

    bool BuildNew(const std::string &InputFile, bool StarFromBare);
    bool Load();
    void Save(); //Save everything in EnvMonteCarlo
    void Reset(ParameterMap);

  private:
    std::string _ConfigFile();
};

class EnvDyson : public Environment {
  public:
    EnvDyson(int pid);
    bool BuildNew(const std::string &InputFile);

  private:
    std::string _FinalQuanFile();
    std::string _FinalStatisFile();
};

int TestEnvironment();
#endif /* defined(__Feynman_Simulator__environment__) */
