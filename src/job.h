//
//  parameter.h
//  Fermion_Simulator
//

#ifndef __Fermion_Simulator__parameter__
#define __Fermion_Simulator__parameter__

#include "definition.h"

enum JobType { MC, DYSON };
// parameters needed to define a job
class Jobs {
  public:
    real T;
    int Order;
    int Seed;
    JobType Type;
    void Read();
};

#endif /* defined(__Fermion_Simulator__parameter__) */
