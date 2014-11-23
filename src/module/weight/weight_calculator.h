//
//  weight_calculator.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/22/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__weight_calculator__
#define __Feynman_Simulator__weight_calculator__

#include "weight_matrix.h"

namespace weight {
class G;
class Sigma;
class GCalculator {
  public:
    GCalculator(G &G_);
    void Dyson(Sigma &Sigma_);
    void FFT(FFT_Mode, fft::Dir);

  private:
    G &_G;
};
}

#endif /* defined(__Feynman_Simulator__weight_calculator__) */
