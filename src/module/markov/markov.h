//
//  markov.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/17/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__markov__
#define __Feynman_Simulator__markov__

#include "utility/convention.h"

namespace diag {
class WormClass;
class Diagram;
}
namespace weight {
class Weight;
class G;
class W;
class Sigma;
class Polar;
class Worm;
}
namespace para {
class ParameterMC;
}
class Lattice;
class RandomFactory;

namespace mc {
class Markov {
  public:
    real Beta;
    Lattice *Lat;
    real *OrderWeight;
    diag::WormClass *Worm;
    weight::Sigma *Sigma;
    weight::Polar *Polar;
    weight::G *G;
    weight::W *W;
    weight::Worm *WormWeight;
    RandomFactory *RNG;

    const static int NUpdates = 4;
    real ProbofCall[NUpdates];
    real SumofProbofCall[NUpdates] = {0.0};

    diag::Diagram *Diag;

    bool BuildNew(para::ParameterMC &, diag::Diagram &, weight::Weight &);
    void ReWeight(para::ParameterMC &);
    void Hop(int);
    void CreateWorm();
    void DeleteWorm();
    void MoveWormOnG();
    void MoveWormOnW();

  private:
    int RandomPickDeltaSpin();
    int RandomPickK();
    int RandomPickDir();
};

int TestMarkov();
}
#endif /* defined(__Feynman_Simulator__markov__) */
