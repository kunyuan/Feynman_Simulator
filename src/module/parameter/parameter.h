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
    real Beta;
    uint MaxTauBin;
    int Order;
    int NSublat;

    //derived
    real T;
    Lattice Lat;

    Message GenerateMessage();
    void UpdateWithMessage(const Message&);

protected:
    bool _FromDict(const Dictionary&);
    Dictionary _ToDict();

    Vec<int> L;
};

class ParaMC : public Parameter {
public:
    long long Counter;
    int Toss;
    int Sweep;
    int Seed;
    RandomFactory RNG;
    real WormSpaceReweight;
    real PolarReweight;
    std::vector<real> OrderReWeight;
    std::vector<real> OrderTimeRatio;

    int PrinterTimer;
    int DiskWriterTimer;
    int MessageTimer;
    int ReweightTimer;

    bool BuildNew(const std::string& InputFile);
    bool FromDict(const Dictionary&);
    Dictionary ToDict();
    void SetTest();
};
}
#endif /* defined(__Feynman_Simulator__state__) */
