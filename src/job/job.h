//
//  parameter.h
//  Fermion_Simulator
//

#ifndef __Fermion_Simulator__parameter__
#define __Fermion_Simulator__parameter__

#include <string>
#include "../utility/convention.h"
#include "../lattice/vector.h"

// parameters needed to define a job
class Jobs {
  private:
    real _T;

  public:
    enum JobType { MC,
                   DYSON };
    int PID;
    Vec<int> L;
    real InitialBeta;
    real DeltaBeta;
    real Beta;
    int Order;
    bool DoesLoad;
    std::string StateFile;
    JobType Type;
    void Read();
};

class JobsMC : public Jobs {
    int Toss;
    int Sample;
    int Sweep;
    int Seed;
    real WormSpaceReWeight;
    real OrderReweight[];
};

class JobsDyson : public Jobs {
};
#endif /* defined(__Fermion_Simulator__parameter__) */
