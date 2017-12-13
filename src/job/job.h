//
//  job.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/1/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__job__
#define __Feynman_Simulator__job__

#include "utility/convention.h"
#include <string>
#include <set>

namespace para {
class Job {
public:
    typedef std::string type;
    std::set<std::string> TypeName = { "MC", "DiagCount" };

    Job(std::string inputfile);
    Job(type, bool, bool, int);

    type Type;
    bool DoesLoad;
    int Sample;
    int PID;
    std::string WeightFile;
    std::string MessageFile;
    std::string StatisticsFile;
    std::string ParaFile;
    std::string LogFile;
    std::string InputFile;
};
}
#endif /* defined(__Feynman_Simulator__job__) */
