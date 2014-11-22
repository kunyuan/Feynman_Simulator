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

class G : public Basic {
  public:
    G(const Lattice &lat, real beta, int order,
      const std::vector<real> &hopping = {0},
      const std::vector<real> &RealChemicalPotential = {0.0, 0.0},
      real ExternalField = 0.0, bool IsTauSymmetric = false);
    Array::array4<Complex> BareWeight;
    Array::array3<Complex> MeasureWeight;
    //Monte Carlo interface
    Complex Weight(const Site &, const Site &, real, real, spin, spin, bool);
    Complex Weight(int, const Site &, const Site &, real, real, spin, spin, bool);
    void Initial(model);
    //Dyson interface

  protected:
    void _InitialTest();
    void _InitialDiagCounter();
    void _InitialBareSpin();
    void _InitialBareHubbardSquare();
    std::vector<real> _Hopping;
    std::vector<real> _RealChemicalPotential;
    real _ExternalField;
};
}

#endif /* defined(__Feynman_Simulator__component__) */
