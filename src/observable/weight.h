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
    Base(const Lattice &, real Beta, int Order);
    std::string _Name;
    real _Beta;
    Lattice _Lat;
    int _Order;
    inline int SpinIndex(spin SpinIn, spin SpinOut);
    inline int SpinIndex(spin *TwoSpinIn, spin *TwoSpinOut);
    inline int TauToBin(real tau);
    inline real BinToTau(int Bin);
    const int MAX_BIN = 128;
    enum Dim {
        ORDER,
        SP,
        SUB,
        VOL,
        TAU
    };
};
/**
*  Estimate gives a estimation of certain quantity with it's error bar'
*/

//TODO: Add fitting function here
class Sigma : Base {
  private:
    real _Norm;
    unsigned int _Shape[5];
    Array::array4<Complex> *_Weight;
    Array::array5<Complex> *_WeightAccu;
    Array::array2<Complex> *_WeightSquareAccu;

  public:
    Sigma(const Lattice &, real Beta, int order);
    ~Sigma();
    int OrderAcceptable(int StartFromOrder, real ErrorThreshold);
    void UpdateWeight(int UpToOrder);
    Complex Weight(const Distance &dR, real dtau, spin, spin);
    Complex WeightOfDelta(spin, spin);
    //TODO: no r dependence, really?
    //    Estimate<Complex> WeightWithError(const Distance &dR, real dtau, spin, spin);
    void Measure(const Distance &, real dtau, spin, spin, int Order, const Complex &);
    void ClearStatistics();
    void SqueezeStatistics(real factor);
    //    std::string PrettyString();
    void SaveState(const std::string &FileName, const std::string &Mode = "a");
    bool LoadState(const std::string &);
};

class Polar : Base {
  private:
    real _Norm;
    unsigned int _Shape[5];
    Array::array4<Complex> *_Weight;
    Array::array5<Complex> *_WeightAccu;
    Array::array2<Complex> *_WeightSquareAccu;

  public:
    Polar(const Lattice &, real Beta, int order);
    ~Polar();
    inline Complex Weight(const Distance &dR, real dtau, spin *, spin *);
    Estimate<Complex> WeightWithError(const Distance &dR, real dtau, spin *, spin *);
};

class G : Base {
  private:
    Array::array4<Complex> *_Weight;
    
  public:
    G(const Lattice &, real Beta, int order);
    ~G();
    void UpdateWeight(int UpToOrder);
    Complex Weight(const Distance &dR, real dtau, spin, spin);
    Complex BareWeight(const Distance &dR, real dtau, spin, spin);
    void SaveState(const std::string &FileName, const std::string &Mode = "a");
    bool LoadState(const std::string &);
};

class W : Base {
  private:
    Array::array4<Complex> *_Weight;

  public:
    W(const Lattice &, real Beta, int order);
    ~W();
    void UpdateWeight(int UpToOrder);
    Complex Weight(const Distance &dR, real dtau, spin *, spin *, bool);
    Complex WeightOfDelta(const Distance &dR, spin *, spin *, bool);
    Complex BareWeight(const Distance &dR, real dtau, spin *, spin *, bool);
    void SaveState(const std::string &FileName, const std::string &Mode = "a");
    bool LoadState(const std::string &);
};
    
class Worm
{
  public:
    real Weight(const Distance &dR, real dtau);
};
    
}

int TestObservable();
#endif /* defined(__Feynman_Simulator__observable__) */
