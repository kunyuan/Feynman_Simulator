//
//  weight_MonteCarlo.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/23/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__weight_MonteCarlo__
#define __Feynman_Simulator__weight_MonteCarlo__

#include "module/weight/component.h"

namespace mc {
namespace weight {
class G : public weight0::G {
  public:
    G(model Model, const Lattice &lat, real beta,
      const std::vector<real> &hopping = {0},
      const std::vector<real> &RealChemicalPotential = {0.0, 0.0},
      real ExternalField = 0.0, weight0::TauSymmetry TauSymmetry = weight0::TauAntiSymmetric);

    Complex Weight(const Site &, const Site &, real, real, spin, spin, bool);
    Complex Weight(int, const Site &, const Site &, real, real, spin, spin, bool);

  private:
    weight0::SmoothTMatrix _MeasureWeight;
};
}
}

#endif /* defined(__Feynman_Simulator__weight_MonteCarlo__) */
