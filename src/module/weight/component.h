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
#include "weight_matrix.h"
namespace weight {
class Sigma;
class GCalculator;
class GMonteCarlo;

class G : public Basic {
  public:
    G(model Model_, const Lattice &lat, real beta,
      const std::vector<real> &hopping = {0},
      const std::vector<real> &RealChemicalPotential = {0.0, 0.0},
      real ExternalField = 0.0, bool IsTauSymmetric = false);

    void BuildNew();
    bool Load(const std::string &FileName);
    void Save(const std::string &FileName, const std::string Mode = "a");
    GCalculator getCalculator();
    GMonteCarlo getMonteCarlo();

  private:
    std::vector<real> _Hopping;
    std::vector<real> _RealChemicalPotential;
    real _ExternalField;

    SmoothTMatrix _BareWeight;
    SmoothTMatrix _SmoothTWeight;
    SmoothTMatrix _MeasureWeight;

    friend class GInitializer;
    friend class GCalculator;
    friend class GMonteCarlo;
};
}

#endif /* defined(__Feynman_Simulator__component__) */
