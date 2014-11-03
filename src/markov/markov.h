//
//  markov.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/17/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__markov__
#define __Feynman_Simulator__markov__

#include "../observable/weight.h"
#include "../parameter/parameter.h"
#include "../diagram/diagram.h"

class Markov {
  public:
    real Beta;
    Lattice *Lat;
    real *OrderWeight;
    WormClass *Worm;
    weight::Sigma *Sigma;
    weight::Polar *Polar;
    weight::G *G;
    weight::W *W;
    weight::Worm *WormWeight;

    const static int NUpdates = 4;
    real ProbofCall[NUpdates];
    real SumofProbofCall[NUpdates]={0.0};

    Diagram *Diag;

    bool BuildNew(ParameterMC &, Diagram &, weight::Weight &);
    void Hop(int);
    void CreateWorm();
    void DeleteWorm();
    void MoveWormOnG();
    void MoveWormOnW();
};

int TestMarkov();
#endif /* defined(__Feynman_Simulator__markov__) */
