//
//  environment.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/17/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__environment__
#define __Feynman_Simulator__environment__

#include <list>
#include "../diagram/diagram.h"
#include "../utility/rng.h"
#include "../observable/weight.h"
#include "../lattice/lattice.h"
#include "../utility/scopeguard.h"

enum JobType { MC,
               DYSON };

JobType GetJobsType(std::string);

class Environment {
  public:
    Environment();
    ~Environment();
    //can only be read from InputFile
    JobType Type;
    bool DoesLoad;
    bool StartFromBare;
    int PID;

    //can be read from StateFile or InputFile
    real Jcp;
    real InitialBeta;
    real DeltaBeta;
    real FinalBeta;
    int Order;

    //derived
    real Beta;
    real T;
    Lattice *Lat;
    Weight::Sigma *Sigma;
    Weight::Polar *Polar;
    Weight::W *W;
    Weight::G *G;

    bool BuildFromFile(std::string InputFile);

    bool LoadState();
    void SaveState(std::string Mode = "a");

    bool LoadGWweight();                       //Load the weight of G,W
    void SaveGWweight(std::string Mode = "a"); //Save the weight of G,W

  protected:
    std::list<std::string> _para;
    Vec<int> _L;
    bool _ReadToPara(std::string);
    std::string _StateFile();
    std::string _ControlFile();
    std::string _GWweightFile();
    std::string _LogFile();
};

class EnvMonteCarlo : public Environment {
  public:
    EnvMonteCarlo();
    ~EnvMonteCarlo();
    //can be read from StateFile or InputFile
    int Version;
    long long Counter;
    int Toss;
    int Sample;
    int Sweep;
    int Seed;
    RandomFactory RNG;
    real WormSpaceReweight;

    Diagram Diag;
    real OrderWeight[MAX_ORDER];

    EstimatorBundle<Complex> cEstimator;
    EstimatorBundle<real> rEstimator;
    EstimatorBundle<real> DetailBalanceEstimator;

    Estimator<real> ZeroOrderWeight;
    Weight::Worm WormWeight;

    bool BuildFromFile(std::string InputFile);
    void Reset();

    bool Load(); //Load according to the flag DoesLoad
    void Save(); //Save everything in EnvMonteCarlo

    bool LoadState();
    void SaveState(std::string Mode = "a");

    bool LoadConfig();
    void SaveConfig(std::string Mode = "a");

    bool LoadStatis();
    void SaveStatis(std::string Mode = "a");

  private:
    void _ReadOrderWeight(char sep = ',');
    void _WriteOrderWeight(ostream &, char sep = ',');
    std::string _ConfigFile();
    std::string _StatisFile();
};

class EnvDyson : public Environment {
  public:
    EnvDyson();
    bool BuildFromFile(std::string InputFile);

  private:
    std::string _TotalStatisFile();
    std::string _FinalQuanFile();
    std::string _FinalStatisFile();
};

#ifndef GET
#define GET(para, thing)                             \
    {                                                \
        stringstream ss(para.front());               \
        (ss) >> (thing);                             \
        if ((ss).fail())                             \
            ABORT("Fail to read " << #thing << "!"); \
        para.pop_front();                            \
    }
#endif

#ifndef PUT
#define PUT(os, thing)                                     \
    {                                                      \
        (os) << (thing) << "    #" << #thing << std::endl; \
        if ((os).fail())                                   \
            ABORT("Fail to write " << #thing << "!");      \
    }
#endif

int TestEnvironment();
#endif /* defined(__Feynman_Simulator__environment__) */
