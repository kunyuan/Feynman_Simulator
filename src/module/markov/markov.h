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
class GClass;
class WClass;
class ExtPoint;
class GammaGClass;
class GammaWClass;
class SigmaClass;
class PolarClass;
class Norm;
}
namespace para {
class ParaMC;
}
class Lattice;
class Site;
class RandomFactory;
class Momentum;

namespace mc {
const int MCUpdates = 19;
const int Gamma3Updates = 4;
const int NUpdates = 23;
class Markov {
public:
    long long* Counter;
    real Beta;
    int Order;
    Lattice* Lat;
    bool runGamma3;
    real* OrderReWeight;
    real* WormSpaceReweight;
    real* PolarReweight;
    real* GammaGReweight;
    real* GammaWReweight;
    diag::Diagram* Diag;
    diag::WormClass* Worm;
    weight::ExtPoint* UExt;
    weight::SigmaClass* Sigma;
    weight::PolarClass* Polar;
    weight::GammaGClass* GammaG;
    weight::GammaWClass* GammaW;
    weight::GClass* G;
    weight::WClass* W;
    RandomFactory* RNG;

    bool BuildNew(para::ParaMC&, diag::Diagram&, weight::Weight&);
    void Reset(para::ParaMC&, diag::Diagram&, weight::Weight&);
    void Hop(int);
    void PrintDetailBalanceInfo();

    void CreateWorm();
    void DeleteWorm();
    void MoveWormOnG();
    void MoveWormOnW();
    void Reconnect();
    void AddInteraction();
    void DeleteInteraction();
    void AddDeltaInteraction();
    void DeleteDeltaInteraction();
    void JumpToOrder0();
    void JumpBackToOrder1();
    void ChangeTauOnVertex();
    void ChangeROnVertex();
    void ChangeRLoop();
    void ChangeMeasureFromGToW();
    void ChangeMeasureFromWToG();
    void ChangeDeltaToContinuous();
    void ChangeContinuousToDelta();
    void ChangeSpinOnVertex();

    //Extra updates for Gamma3
    void JumpToGGGammaG();
    void JumpFromGGGammaGToG();
    void JumpToWWGammaW();
    void JumpFromWWGammaWToW();

private:
    real ProbofCall[NUpdates];
    real SumofProbofCall[NUpdates];
    std::string OperationName[NUpdates];
    real Accepted[NUpdates][MAX_ORDER];
    real Proposed[NUpdates][MAX_ORDER];

    int RandomPickDeltaSpin();
    spin RandomPickSpin();
    Momentum RandomPickK();
    int RandomPickDir();
    real RandomPickTau();
    real ProbTau(real);
    Site RandomPickSite();
    Site RandomPickNeighborSite(const Site&);
    real ProbSite(const Site&);
    real ProbNeighborSite(const Site&);
    bool isNeighbor(const Site&, const Site&);
    bool RandomPickBool();
    enum Operations {
        CREATE_WORM = 0,
        DELETE_WORM,
        MOVE_WORM_G,
        MOVE_WORM_W,
        RECONNECT,
        ADD_INTERACTION,
        DEL_INTERACTION,
        ADD_DELTA_INTERACTION,
        DEL_DELTA_INTERACTION,
        CHANGE_TAU_VERTEX,
        CHANGE_R_VERTEX,
        CHANGE_R_LOOP,
        CHANGE_MEASURE_G2W,
        CHANGE_MEASURE_W2G,
        CHANGE_DELTA2CONTINUS,
        CHANGE_CONTINUS2DELTA,
        CHANGE_SPIN_VERTEX,
        JUMP_TO_ORDER0,
        JUMP_BACK_TO_ORDER1,
        //Extra updates for Gamma3
        JUMP_TO_GGGAMMAG,
        JUMP_FROM_GGGAMMAG_TO_G,
        JUMP_TO_WWGAMMAW,
        JUMP_FROM_WWGAMMAW_TO_W,
        END
    };
    std::string _DetailBalanceStr(Operations op);
    std::string _CheckBalance(Operations op1, Operations op2);
    void _Initial(para::ParaMC&, diag::Diagram&, weight::Weight&);
};

int TestMarkov();
int TestDiagCounter();
}
#endif /* defined(__Feynman_Simulator__markov__) */
