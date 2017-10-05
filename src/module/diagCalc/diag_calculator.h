//
// Created by Kun Chen on 04/10/2017.
//

#ifndef FEYNMANSIMULATOR_DIAG_CALCULATOR_H
#define FEYNMANSIMULATOR_DIAG_CALCULATOR_H

#include "utility/rng.h"
#include "utility/convention.h"
#include "utility/complex.h"
#include "module/weight/weight.h"
#include "estimator/estimator.h"
#include <array>
#include <vector>

namespace weight {
    class Weight;
    class GClass;
    class WClass;
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

#define DiagramsList std::vector<std::array<int, 2*MAX_ORDER>>
#define SpinsList std::vector<std::array<int, 2*MAX_ORDER>>
#define FermiSignList std::vector<int>

namespace diagCalc {

    class DiagramDict{
    public:
        DiagramDict();
        bool FromDict(const Dictionary&);

        std::vector<DiagramsList> AllDiagramConfig;
        std::vector<SpinsList> AllSpinConfig;
        std::vector<FermiSignList> AllFermiSignConfig;
    };

    const int NUpdates = 8;

    class MonteCarloCalc{
    public:
        MonteCarloCalc();

        long long* Counter;
        Lattice* Lat;
        real* OrderReWeight;
        real* PolarReweight;
        diagCalc::DiagramDict *Diag;

        weight::SigmaClass* Sigma;
        weight::PolarClass* Polar;
        weight::GClass* G;
        weight::WClass* W;

        std::array<std::array<Complex, 2*MAX_ORDER>, 2*MAX_ORDER> G_Weight;
        std::array<std::array<Complex, MAX_ORDER>, MAX_ORDER> W_Weight;


        bool BuildNew(para::ParaMC&, diagCalc::DiagramDict &, weight::Weight&);
        void Reset(para::ParaMC&, diagCalc::DiagramDict &, weight::Weight&);
        void Hop(int);
        void CalculateGWTable();
        void PrintDetailBalanceInfo();


        void ChangeMeasureFromGToW();
        void ChangeMeasureFromWToG();
        void ChangeDeltaToContinuous();
        void ChangeContinuousToDelta();

    private:
        real ProbofCall[NUpdates];
        real SumofProbofCall[NUpdates];
        std::string OperationName[NUpdates];
        real Accepted[NUpdates][MAX_ORDER];
        real Proposed[NUpdates][MAX_ORDER];

        spin RandomPickSpin();
        real RandomPickTau();
        real ProbTau(real);
        Site RandomPickSite();
        real ProbSite(const Site&);
        enum Operations {
            CHANGE_R = 0,
            CHANGE_TAU,
            CHANGE_DELTA2CONTINUS,
            CHANGE_CONTINUS2DELTA,
            CHANGE_MEASURE_G2W,
            CHANGE_MEASURE_W2G,
            ADD_INTERACTION,
            DEL_INTERACTION,
            END
        };
        std::string _DetailBalanceStr(Operations op);
        std::string _CheckBalance(Operations op1, Operations op2);
        void _Initial(para::ParaMC&, diagCalc::DiagramDict &, weight::Weight&);
    };

    class MarkovMonitor {
    public:
        MarkovMonitor();

        para::ParaMC *Para;
        diagCalc::DiagramDict *Diag;
        weight::Weight *Weight;

        EstimatorBundle<real> WormEstimator;
        EstimatorBundle<real> PhyEstimator;
        Estimator<real> SigmaEstimator, PolarEstimator;

        bool BuildNew(para::ParaMC &, diagCalc::DiagramDict &, weight::Weight &);
        bool FromDict(const Dictionary &, para::ParaMC &, diagCalc::DiagramDict &, weight::Weight &);
        Dictionary ToDict();

        void Reset(para::ParaMC &, diagCalc::DiagramDict &, weight::Weight &);
        void SqueezeStatistics(real factor);
        bool AdjustOrderReWeight();
        void Measure();
        void AddStatistics();
    };

}

#endif //FEYNMANSIMULATOR_DIAG_CALCULATOR_H
