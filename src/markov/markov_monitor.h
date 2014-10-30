//
//  measure.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/19/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__measure__
#define __Feynman_Simulator__measure__

#include "../environment/environment.h"

class MarkovMonitor {
  private:
    EnvMonteCarlo *_Env;

  public:
    MarkovMonitor(EnvMonteCarlo *);
    void Annealing();
    void SqueezeStatistics();
    void ReWeightEachOrder();
    void Measure();
    void AddStatistics();
};
#endif /* defined(__Feynman_Simulator__measure__) */
