//
//  observable.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/13/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__observable__
#define __Feynman_Simulator__observable__

#include "complex.h"
#include "convention.h"
#include "estimate.h"
#include <iostream>
#include "array.h"

namespace Weight {

class Base {
  protected:
    Base(const Lattice &, real Beta, const std::string &);
    real _Beta;
    Lattice _Lattice;
    std::string _File;
    inline int SpinIndex(spin SpinIn, spin SpinOut);
    inline int SpinIndex(spin *TwoSpinIn, spin *TwoSpinOut);
    inline int TauToBin(real tau);
    inline real BinToTau(int Bin);
    const int MAX_BIN = 128;
};
/**
*  Estimate gives a estimation of certain quantity with it's error bar'
*/

//TODO: Add fitting function here
//TODO: Add reweighting
class Sigma : Base {
  private:
    real _Norm;
    Array::array4<Complex> *_Weight;
    Array::array4<Complex> *_WeightAccu;
    Array::array4<Complex> *_WeightSquareAccu;

  public:
    Sigma(const Lattice &, real Beta, const std::string &);
    ~Sigma();
    void UpdateWeight();
    inline Complex Weight(const Distance &dR, real dtau, spin, spin);
    inline Complex Weight(spin, spin);
    //TODO: no r dependence, really?
    Estimate<Complex> WeightWithError(const Distance &dR, real dtau, spin, spin);
    inline void Measure(const Complex &weight, const Distance &dR, real dtau, spin, spin);
    //    std::string PrettyString();
    void SaveState(std::string);
    bool ReadState(std::string);
};

class Pi : Base {
  private:
    real _Norm;
    Array::array4<Complex> *_Weight;
    Array::array4<Complex> *_WeightAccu;
    Array::array4<Complex> *_WeightSquareAccu;

  public:
    Pi(const Lattice &, real Beta, const std::string &);
    ~Pi();
    void UpdateWeight();
    inline Complex Weight(const Distance &dR, real dtau, spin *, spin *);
    Estimate<Complex> WeightWithError(const Distance &dR, real dtau, spin *, spin *);
    inline void Measure(const Complex &weight, const Distance &dR, real dtau, spin *, spin *);
    //    std::string PrettyString();
    void SaveState(std::string);
    bool ReadState(std::string);
};

class G : Base {
  private:
    Array::array4<Complex> *_Weight;

  public:
    G(const Lattice &, real Beta, const std::string &);
    ~G();
    inline Complex Weight(const Distance &dR, real dtau, spin, spin);
    inline Complex BareWeight(const Distance &dR, real dtau, spin, spin);
    void SaveState(std::string);
    bool ReadState(std::string);
};

class W : Base {
  private:
    Array::array4<Complex> *_Weight;

  public:
    W(const Lattice &, real Beta, const std::string &);
    ~W();
    inline Complex Weight(const Distance &dR, real dtau, spin *, spin *);
    inline Complex Weight(const Distance &dR, spin *, spin *);
    inline Complex BareWeight(const Distance &dR, real dtau, spin *, spin *);
    void SaveState(std::string);
    bool ReadState(std::string);
};
}

int TestObservable();
#endif /* defined(__Feynman_Simulator__observable__) */
