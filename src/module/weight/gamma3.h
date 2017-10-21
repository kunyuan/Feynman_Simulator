//
// Created by kun on 10/18/17.
//

#ifndef FEYNMANSIMULATOR_GAMMA3_H
#define FEYNMANSIMULATOR_GAMMA3_H

#endif //FEYNMANSIMULATOR_GAMMA3_H

#include "index_map.h"
#include "weight_array.h"
#include "weight_estimator.h"
#include "utility/complex.h"
//#include "component.h"

namespace weight {
//     r1,t1  ---\
//                --- 0,0,UP
//  r1,t1+t2  ---/
    class GammaGClass {
        //Only contains the smooth part of the GammaG object
    public:
        GammaGClass(const Lattice &, real Beta, uint MaxTauBin);

//        void BuildNew(GClass *);
//    void BuildTest();

        bool StatisFromDict(const Dictionary &);
        bool WeightFromDict(const Dictionary &);

        Dictionary StatisToDict();
        Dictionary WeightToDict();

        Complex Weight(const Site &, const Site &, const Site &, real, real, real,
                       spin, spin, spin) const;

        void Measure(const Site &, const Site &, const Site &, real, real, real, spin, spin, spin,
                     const Complex &);

        void MeasureNorm(real weight);

        void ClearStatistics();
        void SqueezeStatistics(real factor);

    protected:
        IndexMapSPIN2 _Map;

        WeightArray<4> _WeightAccu; //Gspin, G_r1, Gtau1,dtau2
        WeightArray<4> _Weight; //Gspin, G_r1, Gtau1,dtau2
        uint _CacheIndex[4];
        //don't forget to add template class in weight_array.cpp file
        real _Beta;
        real _Norm; //The normalization factor
        real _NormAccu; //The normalization accumulation
        uint _WeightSize;
        real _dBetaInverse;
        uint _MaxTauBin;
    };

    class GammaWClass {
    public:
        GammaWClass(const Lattice &, real Beta, uint MaxTauBin);

//        void BuildNew();
//    void BuildTest();

        bool StatisFromDict(const Dictionary &);
        bool WeightFromDict(const Dictionary &);
        Dictionary StatisToDict();
        Dictionary WeightToDict();

        Complex Weight(const Site &, const Site &, const Site &, real, real, real,
                       spin*, spin*, spin) const;

        void Measure(const Site &, const Site &, const Site &, real, real, real, spin*, spin*, spin,
                     const Complex &);

        void MeasureNorm(real weight);

        void ClearStatistics();

        void SqueezeStatistics(real factor);

    protected:
        IndexMapSPIN4 _Map;
        WeightArray<5> _Weight; //Wspin, W_r1, W_dr, Wtau1,Wtau2
        WeightArray<5> _WeightAccu; //Wspin, W_r1, W_dr, Wtau1,dtau2
        uint _CacheIndex[5];
        //don't forget to add template class in weight_array.cpp file
        real _Beta;
        real _Norm; //The normalization factor
        real _NormAccu; //The normalization accumulation
        uint _WeightSize;
        real _dBetaInverse;
        uint _MaxTauBin;
        uint _SpinIndex(spin *, spin*) const;
    };
};
