//
//  weight_MC.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/22/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__weight_MC__
#define __Feynman_Simulator__weight_MC__
#include "lattice/lattice.h"
#include "utility/complex.h"

namespace weight {
class G;
class Sigma;
class GMonteCarlo {
  public:
    GMonteCarlo(G &G_);
    Complex Weight(const Site &, const Site &, real, real, spin, spin, bool);
    Complex Weight(int, const Site &, const Site &, real, real, spin, spin, bool);

  private:
    G &_G;
};
}
#endif /* defined(__Feynman_Simulator__weight_MC__) */
