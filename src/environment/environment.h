//
//  environment.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/17/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__environment__
#define __Feynman_Simulator__environment__

#include "module/parameter/parameter.h"
#include "module/diagram/diagram.h"
#include "module/weight/weight.h"
#include "module/markov/markov_monitor.h"
#include "module/markov/markov.h"
#include "job/job.h"

class EnvMonteCarlo {
public:
    EnvMonteCarlo(const para::Job& job, bool IsAllTauSymmetric = false);

    //can be read from StateFile or InputFile
    para::Job Job;
    para::ParaMC Para;
    weight::Weight Weight;
    diag::Diagram Diag;
    mc::Markov Markov;
    mc::MarkovMonitor MarkovMonitor;

    bool BuildNew();
    bool Load();
    void Save(); //Save everything in EnvMonteCarlo
    void DeleteSavedFiles();
    void AdjustOrderReWeight();

    bool ListenToMessage();

private:
    std::string _DiagramFile;
};

int TestEnvironment();
#endif /* defined(__Feynman_Simulator__environment__) */
