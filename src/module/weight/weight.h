//
//  observable.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/13/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__observable__
#define __Feynman_Simulator__observable__

//#include "weight_inherit.h"
#include <string>
#include "utility/convention.h"

class Dictionary;
namespace para {
class ParaMC;
}

namespace weight {

typedef const int flag;
flag SigmaPolar = 1;
flag GW = 2;
class SigmaClass;
class PolarClass;
class GClass;
class WClass;

class Weight {

public:
    Weight(bool IsAllSymmetric = false);
    ~Weight();
    bool _IsAllSymmetric;
    SigmaClass* Sigma;
    PolarClass* Polar;
    GClass* G;
    WClass* W;

    void SetTest(const para::ParaMC&);
    void SetDiagCounter(const para::ParaMC&);
    bool BuildNew(flag, const para::ParaMC&);
    bool FromDict(const Dictionary&, flag, const para::ParaMC&);
    Dictionary ToDict(flag);
    void Anneal(const para::ParaMC&);

private:
    void _AllocateGW(const para::ParaMC&);
    void _AllocateSigmaPolar(const para::ParaMC&);
};
}

#endif /* defined(__Feynman_Simulator__observable__) */
