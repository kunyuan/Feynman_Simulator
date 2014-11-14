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
class ParaMC;
}
class Lattice;
class Site;
class RandomFactory;
class Momentum;

namespace mc {
class Markov {
  public:
    real Beta;
    int Order;
    Lattice *Lat;
    real *OrderWeight;
    diag::WormClass *Worm;
    weight::Sigma *Sigma;
    weight::Polar *Polar;
    weight::G *G;
    weight::W *W;
    weight::Worm *WormWeight;
    RandomFactory *RNG;

    const static int NUpdates = 13;
    real ProbofCall[NUpdates];
    real SumofProbofCall[NUpdates] = {0.0};

    diag::Diagram *Diag;

    bool BuildNew(para::ParaMC &, diag::Diagram &, weight::Weight &);
    void ReWeight(para::ParaMC &);
    void Hop(int);
    
    void CreateWorm();
    void DeleteWorm();
    void MoveWormOnG();
    void MoveWormOnW();
    void Reconnect();
    void AddInteraction();
    void DeleteInteraction();
    void ChangeTauOnVertex();
    void ChangeROnVertex();
    void ChangeRLoop();
    void ChangeMeasureFromGToW();
    void ChangeMeasureFromWToG();
    void ChangeDeltaToContinuous();
    void ChangeContinuousToDelta();
    void ChangeSpinOnVertex();

  private:
    int RandomPickDeltaSpin();
    Momentum RandomPickK();
    int RandomPickDir();
    real RandomPickTau();
    real ProbTau(real);
    Site RandomPickSite();
    real ProbSite(const Site&);
    bool RandomPickBool();
};

int TestMarkov();
int TestDiagCounter();
}
#endif /* defined(__Feynman_Simulator__markov__) */
