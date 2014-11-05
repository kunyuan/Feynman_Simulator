//
//  environment.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/17/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__environment__
#define __Feynman_Simulator__environment__

#include "utility/rng.h"
#include "utility/scopeguard.h"
#include "lattice/lattice.h"
#include "module/parameter/parameter.h"
#include "module/diagram/diagram.h"
#include "module/observable/weight.h"
#include "module/markov/markov_monitor.h"
#include "module/markov/markov.h"
#include "module/dyson/dyson.h"

class Environment {
  public:
    int PID;
    enum flag { Bare,
                OldGW };

  protected:
    std::string _ParameterFile;
    std::string _GWweightFile;
    std::string _WeightFile;
    std::string _StatisticsFile;
    Environment(int pid);
};

class EnvMonteCarlo : public Environment {
  public:
    EnvMonteCarlo(int pid);

    //can be read from StateFile or InputFile
    para::ParaMC Para;
    weight::Weight Weight;
    diag::Diagram Diag;
    mc::Markov Grasshopper;
    mc::MarkovMonitor Scarecrow;

    bool BuildNew(const std::string &InputFile, bool StartFromBare);
    bool Load();
    void Save(); //Save everything in EnvMonteCarlo
    bool ReLoad();
    void DeleteSavedFiles();
    bool ReWeight();
    bool CheckStatus();

  private:
    std::string _DiagramFile;
};

class EnvDyson : public Environment {
  public:
    EnvDyson(int pid);

    para::ParaDyson Para;
    weight::Weight Weight;
    dyson::Dyson Dyson;

    bool BuildNew(const std::string &InputFile, bool StartFromBare);
    bool CanBeLoad();
    bool Load();
    void Save();

  private:
};

int TestEnvironment();
#endif /* defined(__Feynman_Simulator__environment__) */
