//
//  parameter.h
//  Fermion_Simulator
//

#ifndef __Fermion_Simulator__parameter__
#define __Fermion_Simulator__parameter__

#include <string>
#include <iostream>
#include "../utility/convention.h"
#include "../lattice/vector.h"

namespace Jobs {

typedef shared_ptr<std::ifstream> jobstream;

enum JobType { MC,
               DYSON };

JobType GetJobsType(std::string);

class JobsBase {
  protected:
    jobstream ifs;

  public:
    JobsBase(std::string InputFile);
    JobType Type;
    int PID;
    Vec<int> L;
    real Jcp;
    real InitialBeta;
    real DeltaBeta;
    real Beta;
    int Order;
    bool DoesLoad;
    std::string StateFile;
    void TestJob();
};

class JobsMC : public JobsBase {
  public:
    JobsMC();
    JobsMC(std::string InputFile);
    int Toss;
    int Sample;
    int Sweep;
    int Seed;
    real WormSpaceReweight;
    real OrderReweight[];
};

class JobsDyson : public JobsBase {
  public:
    JobsDyson();
    JobsDyson(std::string InputFile);
};
}
#endif /* defined(__Fermion_Simulator__parameter__) */
