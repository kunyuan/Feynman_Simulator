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

class WeightNoMeasure {
  protected:
    WeightNoMeasure(const Lattice &, real Beta, int Order, int SpinVol, std::string);
    ~WeightNoMeasure();
    std::string _Name;
    real _Beta;
    real _dBeta;
    real _dBetaInverse;
    int _Order;
    Lattice _Lat;
    unsigned int _Shape[5];
    Array::array4<Complex> *_Weight;

    int SpinIndex(spin SpinIn, spin SpinOut);
    int SpinIndex(spin *TwoSpinIn, spin *TwoSpinOut);
    int TauToBin(real tau);
    real BinToTau(int Bin);
    const int MAX_BIN = 128;
    enum Dim { ORDER,
               SP,
               SUB,
               VOL,
               TAU };

  public:
    void SaveState(const std::string &FileName, const std::string &Mode = "a");
    bool LoadState(const std::string &);
};

class WeightNeedMeasure : public WeightNoMeasure {
  protected:
    real _Norm;
    Array::array5<Complex> *_WeightAccu;
    EstimatorBundle<Complex> _Average;

  public:
    WeightNeedMeasure(const Lattice &, real Beta, int order, int Spine, std::string);
    ~WeightNeedMeasure();

    Estimate<Complex> WeightWithError(int order);

    int OrderAcceptable(int StartFromOrder, real ErrorThreshold);
    void UpdateWeight(int UpToOrder);

    void AddStatistics();
    void ClearStatistics();
    void SqueezeStatistics(real factor);
    //    std::string PrettyString();
    void SaveState(const std::string &FileName, const std::string &Mode = "a");
    bool LoadState(const std::string &);
};

//TODO: Add fitting function here
class Sigma : public WeightNeedMeasure {
  public:
    Sigma(const Lattice &, real Beta, int order);
    Complex Weight(const Distance &dR, real dtau, spin, spin);
    Complex WeightOfDelta(spin, spin);
    void Measure(const Distance &, real dtau, spin, spin, int Order, const Complex &);
};

class Polar : public WeightNeedMeasure {
  public:
    Polar(const Lattice &, real Beta, int order);
    Complex Weight(const Distance &dR, real dtau, spin *, spin *);
    void Measure(const Distance &, real dtau, spin *, spin *, int Order, const Complex &);
};

class G : public WeightNeedMeasure {
  public:
    G(const Lattice &, real Beta, int order);
    Complex Weight(const Distance &dR, real dtau, spin, spin);
    Complex BareWeight(const Distance &dR, real dtau, spin, spin);
};

class W : public WeightNeedMeasure {
  public:
    W(const Lattice &, real Beta, int order);
    Complex Weight(const Distance &dR, real dtau, spin *, spin *);
    Complex WeightOfDelta(const Distance &dR, spin *, spin *);
    Complex BareWeight(const Distance &dR, real dtau, spin *, spin *);
};
}

int TestObservable();
#endif /* defined(__Feynman_Simulator__observable__) */
