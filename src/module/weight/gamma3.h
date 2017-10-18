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

namespace weight {
//     r1,t1  ---\
//                --- 0,0,UP
//  r1,t1+t2  ---/
    class GammaGClass {
    public:
        GammaGClass(const Lattice &, real Beta, uint MaxTauBin, real Norm);

        void BuildNew();
//    void BuildTest();

        void Reset(real Beta);

        bool FromDict(const Dictionary &);

        Dictionary ToDict();

        Complex Weight(const Site &, const Site &, const Site &, real, real, real,
                       spin, spin, spin) const;

        void Measure(const Site &, const Site &, const Site &, real, real, real, spin, spin, spin,
                     const Complex &);

        void MeasureNorm(real weight);

        void ClearStatistics();
        void SqueezeStatistics(real factor);

    protected:
        IndexMapSPIN2 _Map;

        WeightArray<6> _WeightAccu; //Gspin, Uspin, Gsub, Usub, G_r1, Gtau1,dtau2
        real _Beta;
        real _Norm; //The normalization factor
        real _NormAccu; //The normalization accumulation
        uint _WeightSize;
    };
};
