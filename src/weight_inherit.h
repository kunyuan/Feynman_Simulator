//
//  weight_inheritance.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/7/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef Feynman_Simulator_weight_inheritance_h
#define Feynman_Simulator_weight_inheritance_h
#include "weight_estimator.h"

namespace weight {

//TODO: Add fitting function here
class Sigma : public WeightNeedMeasure {
  public:
    Sigma(const Lattice &, real Beta, int order);
    Complex Weight(const Site &, const Site &, real, real, spin, spin);
    Complex WeightOfDelta(spin, spin);
    void Measure(const Site &, const Site &, real, real, spin, spin, int Order, const Complex &);
    void FFT(fft::Dir, Mode);
};

class Polar : public WeightNeedMeasure {
  public:
    Polar(const Lattice &, real Beta, int order);
    Complex Weight(const Site &, const Site &, real, real, spin *, spin *);
    void Measure(const Site &, const Site &, real, real, spin *, spin *, int Order, const Complex &);
    void FFT(fft::Dir, Mode);
};
class G : public WeightNoMeasure {
  public:
    G(const Lattice &, real Beta, int order);
    Complex Weight(const Site &, const Site &, real, real, spin, spin, bool);
    Complex Weight(int, const Site &, const Site &, real, real, spin, spin, bool);
    void InitialWithBare();
    void FFT(fft::Dir, Mode);
};

class W : public WeightNoMeasure {
  public:
    W(const Lattice &, real Beta, int order);
    Complex Weight(const Site &, const Site &, real, real, spin *, spin *, bool, bool, bool);
    Complex Weight(int, const Site &, const Site &, real, real, spin *, spin *, bool, bool, bool);
    void InitialWithBare();
    void FFT(fft::Dir, Mode);
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

#endif
