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
#include "weight_estimator.h"
#include "utility/complex.h"

namespace weight {
class G : public Basic {
    friend class GInitializer;

  public:
    G(const Lattice &lat, real beta, TauSymmetry TauSymmetry = TauAntiSymmetric);
    void BuildNew(model Model,
                  const std::vector<real> &Hopping = {0.0},
                  const std::vector<Complex> &ChemicalPotential = {0.0, 0.0},
                  real ExternalField = 0.0);
    void BuildTest();
    void Reset(real Beta);

    Complex Weight(const Site &, const Site &, real, real, spin, spin, bool);
    Complex Weight(int, const Site &, const Site &, real, real, spin, spin, bool);

  private:
    std::vector<real> _Hopping;
    std::vector<Complex> _ChemicalPotential;
    real _ExternalField;
    weight::SmoothTMatrix _MeasureWeight;
    IndexMapSPIN2 _Map;
};

/**
*  W is the interaction. An assumption is made here: translational and \emp{MIRROR} symmetry of the lattice (constructed by unit cells) are imposed on interaction.
    The mirror symmetry is only required on the level of the whole lattice, not within a unit cell.
*/
class W : public Basic {
    friend class WInitializer;

  public:
    W(const Lattice &lat, real Beta);
    void BuildNew(model Model,
                  const std::vector<real> &Interaction = {0.0},
                  real ExternalField = 0.0);
    void BuildTest();
    void WriteBareToASCII();
    void Reset(real Beta);

    Complex Weight(const Site &, const Site &, real, real, spin *, spin *, bool, bool, bool);
    Complex Weight(int, const Site &, const Site &, real, real, spin *, spin *, bool, bool, bool);

  protected:
    std::vector<real> _Interaction;
    real _ExternalField;
    weight::SmoothTMatrix _MeasureWeight;
    IndexMapSPIN4 _Map;
};

class Sigma : public Basic {
  public:
    Sigma(const Lattice &, real Beta, int MaxOrder,
          TauSymmetry Symmetry = TauAntiSymmetric);
    void BuildNew(model Model);
    void BuildTest();

    void Reset(real Beta);
    bool Load(const std::string &FileName);
    void Save(const std::string &FileName, const std::string Mode = "a");

    Complex Weight(const Site &, const Site &, real, real, spin, spin);
    void Measure(const Site &, const Site &, real, real, spin, spin,
                 int Order, const Complex &);
    void UpdateWeight(int order);

    WeightEstimator Estimator;

  protected:
    IndexMapSPIN2 _Map;
};

class Polar : public Basic {
  public:
    Polar(const Lattice &, real Beta, int MaxOrder);
    void BuildNew(model Model);
    void BuildTest();

    void Reset(real Beta);
    bool Load(const std::string &FileName);
    void Save(const std::string &FileName, const std::string Mode = "a");

    Complex Weight(const Site &, const Site &, real, real, spin *, spin *);
    void Measure(const Site &, const Site &, real, real, spin *, spin *,
                 int Order, const Complex &);
    void UpdateWeight(int order);

    WeightEstimator Estimator;

  protected:
    IndexMapSPIN4 _Map;
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

int TestWeight();
}

#endif /* defined(__Feynman_Simulator__component__) */
