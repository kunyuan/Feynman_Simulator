//
//  markov.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/17/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__markov__
#define __Feynman_Simulator__markov__

#include "../environment/environment.h"

class Markov {
  private:
    real Beta;
    Lattice *Lat;
    real *OrderWeight;
    WormClass *Worm;
    Weight::Sigma *Sigma;
    Weight::Polar *Polar;
    Weight::G *G;
    Weight::W *W;
    Weight::Worm *WormWeight;

    const static int NUpdates = 2;
    real ProbofCall[NUpdates];

  public:
    Diagram *Diag;
    Markov(EnvMonteCarlo *);
    void Hop(int);
    void CreateWorm();
    void DeleteWorm();
};

int TestMarkov();
#endif /* defined(__Feynman_Simulator__markov__) */
