//
//  measure.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/19/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__measure__
#define __Feynman_Simulator__measure__

#include "environment.h"

class MarkovMonitor {
  private:
    EnvMoneCarlo *_Env;
    Diagram *_Diag;
    RandomFactory *_RNG;

  public:
    MarkovMonitor(EnvMoneCarlo *);
};
#endif /* defined(__Feynman_Simulator__measure__) */
