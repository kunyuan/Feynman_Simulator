//
//  component.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/22/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__component__
#define __Feynman_Simulator__component__

#include "weight_basic.h"
#include "weight_estimator.h"
#include "utility/complex.h"

namespace weight {
enum model { DiagCount = 0,
             Trivial };

class Worm {
public:
    static real Weight(const Site&, const Site&, real, real)
    {
        return 1.0;
    }
};

class Norm {
public:
    static real Weight()
    {
        return 10.0;
    }
};

class G : public Basic {
public:
    G(const Lattice& lat, real beta, uint MaxTauBin,
      TauSymmetry TauSymmetry = TauAntiSymmetric);
    void BuildTest(weight::model);
    void Reset(real Beta);

    Complex Weight(const Site&, const Site&, real, real, spin, spin, bool) const;
    Complex Weight(int, const Site&, const Site&, real, real, spin, spin, bool) const;

private:
    weight::SmoothTMatrix _MeasureWeight;
    IndexMapSPIN2 _Map;
};

/**
*  W is the interaction. An assumption is made here: translational and \emp{MIRROR} symmetry of the lattice (constructed by unit cells) are imposed on interaction.
    The mirror symmetry is only required on the level of the whole lattice, not within a unit cell.
*/
class W : public Basic {
public:
    W(const Lattice& lat, real Beta, uint MaxTauBin);
    void BuildTest(weight::model);
    void WriteBareToASCII();
    void Reset(real Beta);

    Complex Weight(const Site&, const Site&, real, real, spin*, spin*, bool, bool, bool) const;
    Complex Weight(int, const Site&, const Site&, real, real, spin*, spin*, bool, bool, bool) const;

protected:
    weight::SmoothTMatrix _MeasureWeight;
    IndexMapSPIN4 _Map;
};

class Sigma : public Basic {
public:
    Sigma(const Lattice&, real Beta, uint MaxTauBin, int MaxOrder,
          TauSymmetry Symmetry = TauAntiSymmetric, real Norm = Norm::Weight());
    void BuildNew();
    void BuildTest();

    void Reset(real Beta);
    bool FromDict(const Dictionary&);
    Dictionary ToDict();

    Complex Weight(const Site&, const Site&, real, real, spin, spin) const;
    void Measure(const Site&, const Site&, real, real, spin, spin,
                 int Order, const Complex&);
    void UpdateWeight(int order);

    WeightEstimator Estimator;

protected:
    IndexMapSPIN2 _Map;
};

class Polar : public Basic {
public:
    Polar(const Lattice&, real Beta, uint MaxTauBin, int MaxOrder, real Norm = Norm::Weight());
    void BuildNew();
    void BuildTest();

    void Reset(real Beta);
    bool FromDict(const Dictionary&);
    Dictionary ToDict();

    Complex Weight(const Site&, const Site&, real, real, spin*, spin*) const;
    void Measure(const Site&, const Site&, real, real, spin*, spin*,
                 int Order, const Complex&);
    void UpdateWeight(int order);

    WeightEstimator Estimator;

protected:
    IndexMapSPIN4 _Map;
};

int TestWeight();
}

#endif /* defined(__Feynman_Simulator__component__) */
