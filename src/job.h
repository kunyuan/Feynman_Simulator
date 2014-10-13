//
//  parameter.h
//  Fermion_Simulator
//

#ifndef __Fermion_Simulator__parameter__
#define __Fermion_Simulator__parameter__

#include "definition_global.h"

// parameters needed to define a job
class JobPara {
  public:
    real Temp;
    real Beta;
    int Order;
    int Seed;

    void Read();
};

#endif /* defined(__Fermion_Simulator__parameter__) */
