//
//  markov.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/17/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__markov__
#define __Feynman_Simulator__markov__

#include <string>
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
    long long* Counter;
    real Beta;
    int Order;
    Lattice* Lat;
    real* OrderWeight;
    diag::Diagram* Diag;
    diag::WormClass* Worm;
    weight::Sigma* Sigma;
    weight::Polar* Polar;
    weight::G* G;
    weight::W* W;
    RandomFactory* RNG;

    bool BuildNew(para::ParaMC&, diag::Diagram&, weight::Weight&);
    void ReWeight(para::ParaMC&);
    void Hop(int);
    void PrintDetailBalanceInfo();

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
    const static int NUpdates = 15;
    real ProbofCall[NUpdates] = { 0.0 };
    real SumofProbofCall[NUpdates] = { 0.0 };
    std::string OperationName[NUpdates];
    real Accepted[NUpdates][MAX_ORDER] = { 0.0 };
    real Proposed[NUpdates][MAX_ORDER] = { 0.0 };

    int RandomPickDeltaSpin();
    Momentum RandomPickK();
    int RandomPickDir();
    real RandomPickTau();
    real ProbTau(real);
    Site RandomPickSite();
    real ProbSite(const Site&);
    bool RandomPickBool();
    enum Operations {
        CREATE_WORM = 0,
        DELETE_WORM,
        MOVE_WORM_G,
        MOVE_WORM_W,
        RECONNECT,
        ADD_INTERACTION,
        DEL_INTERACTION,
        CHANGE_TAU_VERTEX,
        CHANGE_R_VERTEX,
        CHANGE_R_LOOP,
        CHANGE_MEASURE_G2W,
        CHANGE_MEASURE_W2G,
        CHANGE_DELTA2CONTINUS,
        CHANGE_CONTINUS2DELTA,
        CHANGE_SPIN_VERTEX
    };
    std::string _DetailBalanceStr(Operations op);
    std::string _CheckBalance(Operations op1, Operations op2);
};

int TestMarkov();
int TestDiagCounter();
}
#endif /* defined(__Feynman_Simulator__markov__) */
