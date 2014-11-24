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
#include "utility/array.h"
#include "utility/complex.h"

namespace weight0 {
typedef Array::array4<Complex> SmoothTMatrix;
typedef Array::array3<Complex> DeltaTMatrix;

class G : public BasicWithTwoSpins {
  public:
    G(model Model, const Lattice &lat, real beta,
      const std::vector<real> &hopping = {0},
      const std::vector<real> &RealChemicalPotential = {0.0, 0.0},
      real ExternalField = 0.0, TauSymmetry TauSymmetry = TauAntiSymmetric);

    void BuildNew();
    void BuildTest();
    bool Load(const std::string &FileName);
    void Save(const std::string &FileName, const std::string Mode = "a");

  protected:
    std::vector<real> _Hopping;
    std::vector<real> _RealChemicalPotential;
    real _ExternalField;

    SmoothTMatrix _BareWeight;
    SmoothTMatrix _SmoothTWeight;

    friend class GInitializer;
};
}

#endif /* defined(__Feynman_Simulator__component__) */
