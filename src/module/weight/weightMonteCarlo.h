//
//  weightMonteCarlo.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/9/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__weightMonteCarlo__
#define __Feynman_Simulator__weightMonteCarlo__

#include "weight_class.h"

namespace weight0 {
namespace mc {

class Sigma : protected WeightEstimator {
  public:
    Sigma(const Lattice &, real Beta, int order);
    void Measure(const Site &, const Site &, real, real, spin, spin, int Order, const Complex &);
    void MeasureNorm();
};

class Polar : protected WeightEstimator {
  public:
    Polar(const Lattice &, real Beta, int order);
    void Measure(const Site &, const Site &, real, real, spin *, spin *, int Order, const Complex &);
    void MeasureNorm();
};

class G : protected WeightArray {
  public:
    G(const Lattice &, real Beta);
    Complex Weight(const Site &, const Site &, real, real, spin, spin, bool);
    Complex Weight(int, const Site &, const Site &, real, real, spin, spin, bool);
    void InitialWithBare();
};

class W : protected WeightArray {
  public:
    W(const Lattice &, real Beta);
    Complex Weight(const Site &, const Site &, real, real, spin *, spin *, bool, bool, bool);
    Complex Weight(int Dir, const Site &, const Site &, real, real, spin *, spin *, bool, bool, bool);
    void InitialWithBare();
};

class Worm {
  public:
    static real Weight(const Site &, const Site &, real, real)
    {
        return 1.0;
    }
};

class Norm {
  public:
    static real Weight()
    {
        return 1.0;
    }
};
}
}

#endif /* defined(__Feynman_Simulator__weightMonteCarlo__) */
