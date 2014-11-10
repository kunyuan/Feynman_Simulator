//
//  job.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/1/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__job__
#define __Feynman_Simulator__job__

#include "utility/parser.h"
#include <set>

namespace para {
class Job {
  public:
    typedef std::string type;
    std::set<std::string> TypeName = {"MC", "DYSON"};

    Job(std::string inputfile);
    Job(type, bool, bool, int);

    type Type;
    bool DoesLoad;
    bool StartFromBare;
    int PID;
    std::string InputFile;

  private:
    SimpleParser _Para;
};
}
#endif /* defined(__Feynman_Simulator__job__) */
