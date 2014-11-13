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
#include "module/weight/weight.h"
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
    EnvMonteCarlo(int pid, bool IsAllTauSymmetric = false);

    //can be read from StateFile or InputFile
    para::ParaMC Para;
    weight::Weight Weight;
    diag::Diagram Diag;
    mc::Markov Grasshopper;
    mc::MarkovMonitor Scarecrow;

    bool BuildNew(const std::string &InputFile, bool StartFromBare);
    bool Load();
    void Save(); //Save everything in EnvMonteCarlo
    void DeleteSavedFiles();
    bool ReWeight();

    bool ListenToMessage();

  private:
    std::string _DiagramFile;
};

class EnvDyson : public Environment {
  public:
    EnvDyson(int pid, bool IsAllTauSymmetric = false);

    para::ParaDyson Para;
    weight::Weight Weight;
    dyson::Dyson Dyson;

    bool BuildNew(const std::string &InputFile, bool StartFromBare);
    bool Load();
    void Save();
    void UpdateWeight();
    //Update the weight of Sigma and Polar according to Para.ErrorThreshold and Para.OrderAccepted
    void BroadcastMessage();

  private:
};

int TestEnvironment();
#endif /* defined(__Feynman_Simulator__environment__) */
