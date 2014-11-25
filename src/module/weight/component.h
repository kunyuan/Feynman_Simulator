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
#include "utility/complex.h"

namespace weight0 {
class G : public Basic {
    friend class GInitializer;

  public:
    G(model Model, const Lattice &lat, real beta,
      const std::vector<real> &Hopping = {0},
      const std::vector<real> &RealChemicalPotential = {0.0, 0.0},
      real ExternalField = 0.0, TauSymmetry TauSymmetry = TauAntiSymmetric);

    void BuildNew();
    void BuildTest();
    void Reset(real Beta);

  protected:
    std::vector<real> _Hopping;
    std::vector<real> _RealChemicalPotential;
    real _ExternalField;

    IndexMapSPIN2 _Map;
};

/**
*  W is the interaction. An assumption is made here: translational and \emp{MIRROR} symmetry of the lattice (constructed by unit cells) are imposed on interaction.
    The mirror symmetry is only required on the level of the whole lattice, not within a unit cell.
*/
class W : public Basic {
    friend class WInitializer;

  public:
    W(model Model, const Lattice &lat, real Beta,
      const std::vector<real> &Interaction, real ExternalField = 0.0);
    void BuildNew();
    void BuildTest();
    void WriteBareToASCII();
    void Reset(real Beta);

  protected:
    std::vector<real> _Interaction;
    real _ExternalField;

    IndexMapSPIN4 _Map;
};

class Sigma : public Basic {
  public:
    Sigma(model Model, const Lattice &, real Beta, TauSymmetry Symmetry = TauAntiSymmetric);
    void Reset(real Beta);

  protected:
    IndexMapSPIN2 _Map;
};

class Polar : public Basic {
  public:
    Polar(model Model, const Lattice &, real Beta);
    void Reset(real Beta);

  protected:
    IndexMapSPIN4 _Map;
};

int TestWeight();
}

#endif /* defined(__Feynman_Simulator__component__) */
