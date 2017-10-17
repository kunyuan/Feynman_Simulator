//
// Created by Kun Chen on 04/10/2017.
//

#ifndef FEYNMANSIMULATOR_DIAG_CALCULATOR_H
#define FEYNMANSIMULATOR_DIAG_CALCULATOR_H

#include "utility/rng.h"
#include "utility/convention.h"
#include "utility/complex.h"
#include "lattice/lattice.h"
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

    class Conf{
    public:
        void Initialize();

        int Order;
        bool MeasureGLine;
        Complex Phase;
        real Weight;

        std::array<Site, 2*MAX_ORDER> R_list;
        std::array<real, 2*MAX_ORDER> Tau_list;
    };

<<<<<<< HEAD
    const int NUpdates = 10;

    class MonteCarlo{
    public:
=======
    class MonteCarlo{
    public:
        MonteCarlo();
>>>>>>> calc

        long long* Counter;
        Lattice* Lat;
        int Order;
        real Beta;
        real* OrderReWeight;
        real* PolarReweight;
        diagCalc::DiagramDict *Diag;

        weight::SigmaClass* Sigma;
        weight::PolarClass* Polar;
        weight::GClass* G;
        weight::WClass* W;
        RandomFactory* RNG;

        std::array<std::array<std::array<Complex, 2*MAX_ORDER>, 2*MAX_ORDER>,2> G_Weight;
        std::array<std::array<std::array<Complex, MAX_ORDER>, MAX_ORDER>, 16> W_Weight;

        std::vector<spin> IndexToSpin;

        Conf DiagConf;

        bool BuildNew(para::ParaMC&, diagCalc::DiagramDict &, weight::Weight&);
        void Reset(para::ParaMC&, diagCalc::DiagramDict &, weight::Weight&);
        void Hop(int);

        void CalculateGWTable();
        Complex SumAllDiagrams();

        void PrintDetailBalanceInfo();

        void ChangeR();
        void ChangeTau();
        void ChangeDeltaToContinuous();
        void ChangeContinuousToDelta();
        void ChangeMeasureFromGToW();
        void ChangeMeasureFromWToG();
        void AddInteraction();
        void DeleteInteraction();
        void JumpToOrder0();
        void JumpToOrder1();

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
            JUMP_TO_ORDER0,
            JUMP_TO_ORDER1,
            END
        };
        std::string _DetailBalanceStr(Operations op);
        std::string _CheckBalance(Operations op1, Operations op2);
        void _Initial(para::ParaMC&, diagCalc::DiagramDict &, weight::Weight&);


        static int SpinIndex(spin SpinInIn, spin SpinInOut, spin SpinOutIn, spin SpinOutOut);
        static int SpinIndex(const spin* TwoSpinIn, const spin* TwoSpinOut);
    };

    class MarkovMonitor {
    public:

        para::ParaMC *Para;
        diagCalc::DiagramDict *Diag;
        diagCalc::Conf *DiagConf;
        weight::Weight *Weight;

//        EstimatorBundle<real> WormEstimator;
        EstimatorBundle<real> PhyEstimator;
        Estimator<real> SigmaEstimator, PolarEstimator;

        bool BuildNew(para::ParaMC &, diagCalc::DiagramDict &, weight::Weight &, diagCalc::Conf &);
        bool FromDict(const Dictionary &, para::ParaMC &, diagCalc::DiagramDict &, weight::Weight &, diagCalc::Conf &);
        Dictionary ToDict();

        void Reset(para::ParaMC &, diagCalc::DiagramDict &, weight::Weight &, diagCalc::Conf &);
        void SqueezeStatistics(real factor);
        bool AdjustOrderReWeight();
        void Measure();
        void AddStatistics();
    };

}

#endif //FEYNMANSIMULATOR_DIAG_CALCULATOR_H
