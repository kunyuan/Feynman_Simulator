//
//  weightDyson.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/9/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__weightDyson__
#define __Feynman_Simulator__weightDyson__

#include "weight_class.h"
namespace weight0 {
namespace dyson {
class Sigma : public WeightArray {
  public:
    Sigma(const Lattice &, real Beta);
    void FFT(fft::Dir, Mode);
};

class Polar : public WeightArray {
  public:
    Polar(const Lattice &, real Beta);
    void FFT(fft::Dir, Mode);
};
class G : public WeightArray {
  public:
    G(const Lattice &, real Beta);
    void FFT(fft::Dir, Mode);
};

class W : public WeightArray {
  public:
    W(const Lattice &, real Beta);
    void FFT(fft::Dir, Mode);
};
}
}

#endif /* defined(__Feynman_Simulator__weightDyson__) */
