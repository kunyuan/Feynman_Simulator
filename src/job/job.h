//
//  parameter.h
//  Fermion_Simulator
//

#ifndef __Fermion_Simulator__parameter__
#define __Fermion_Simulator__parameter__

#include "../utility/convention.h"

// parameters needed to define a job
class Jobs {
  private:
    real _T;

  public:
    real Beta;
    int Order;
    int Seed;
    JobType Type;
    void Read();
};

#endif /* defined(__Fermion_Simulator__parameter__) */
