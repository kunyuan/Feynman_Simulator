//
//  observable.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/13/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__observable__
#define __Feynman_Simulator__observable__

#include "weight_inherit.h"

namespace para {
class Parameter;
}

namespace weight {

typedef const int flag;
flag SigmaPolar = 1;
flag GW = 2;

class Weight {
  public:
    Weight(bool IsAllSymmetric=false);
    ~Weight();
    bool _IsAllSymmetric;
    Sigma *Sigma;
    Polar *Polar;
    G *G;
    W *W;
    Worm WormWeight;

    void SetTest(const para::Parameter &);
    void SetDiagCounter(const para::Parameter &);
    bool BuildNew(flag, const para::Parameter &);
    bool Load(const std::string &InputFile, flag, const para::Parameter &);
    void Save(const std::string &InputFile, flag, string Mode = "a");
    void ReWeight(flag, const para::Parameter &);
    int UpdateSigmaPolarWeight(int OrderAccepted, real ErrorThreshold);

  private:
    void _AllocateGW(const para::Parameter &);
    void _AllocateSigmaPolar(const para::Parameter &);
};
int TestWeight();
}

#endif /* defined(__Feynman_Simulator__observable__) */
