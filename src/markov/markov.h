//
//  markov.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/17/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__markov__
#define __Feynman_Simulator__markov__

#include "environment.h"
#include "weight.h"
#include "rng.h"

class Markov {
  private:
    real Beta;
    Lattice *Lat;
    real *OrderWeight;
    Diagram *Diag;
    Weight::Sigma *Sigma;
    Weight::Polar *Polar;

  public:
    Markov(EnvMonteCarlo *);
    void Hop(int &&);
    void CreateWorm();
    void DeleteWorm();
};

#endif /* defined(__Feynman_Simulator__markov__) */
