//
//  state.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/2/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__state__
#define __Feynman_Simulator__state__

#include "lattice/lattice.h"
#include "message.h"
#include "utility/parser.h"
#include "utility/rng.h"
#include "utility/convention.h"
#include <string>

namespace para {

class Parameter {
  public:
    int Version;
    std::vector<real> Interaction;
    real ExternalField;
    real InitialBeta;
    real DeltaBeta;
    real FinalBeta;
    int Order;

    //derived
    real Beta;
    real T;
    Lattice Lat;

    Message GenerateMessage();
    void UpdateWithMessage(const Message &);

  protected:
    bool _BuildNew(const std::string &InputFile);
    bool _Load(const std::string &InputFile);
    void _SavePreparation();

    SimpleParser _para;
    Vec<int> L;
};

class ParaMC : public Parameter {
  public:
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
    void SetTest();
};

class ParaDyson : public Parameter {
  public:
    int OrderAccepted;
    real ErrorThreshold;
    int SleepTime;

    bool BuildNew(const std::string &InputFile);
    bool Load(const std::string &InputFile);
    void Save(const std::string &InputFile, std::string Mode = "a");
    void SetTest();
};
}
#endif /* defined(__Feynman_Simulator__state__) */
