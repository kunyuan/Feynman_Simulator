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
#include "utility/rng.h"
#include "utility/complex.h"
#include "utility/convention.h"
#include <string>

class Dictionary;
namespace para {

class Parameter {
public:
    int Version;
    real InitialBeta;
    real DeltaBeta;
    real FinalBeta;
    int Order;
    int NSublat;

    //derived
    real Beta;
    real T;
    Lattice Lat;

    Message GenerateMessage();
    void UpdateWithMessage(const Message&);

protected:
    bool _BuildNew(const std::string& InputFile);
    bool _Load(const std::string& InputFile);
    bool _FromDict(const Dictionary&);
    Dictionary _ToDict();

    Vec<int> L;
};

class ParaMC : public Parameter {
public:
    long long Counter;
    uint MaxTauBin;
    int Toss;
    int Sample;
    int Sweep;
    int Seed;
    RandomFactory RNG;
    real WormSpaceReweight;
    std::vector<real> OrderReWeight;

    bool BuildNew(const std::string& InputFile);
    bool Load(const std::string& InputFile);
    void Save(const std::string& InputFile, std::string Mode = "a");
    bool FromDict(const Dictionary&);
    Dictionary ToDict();
    void SetTest();
};
}
#endif /* defined(__Feynman_Simulator__state__) */
