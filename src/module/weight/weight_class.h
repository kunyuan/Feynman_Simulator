//
//  weight_array.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/8/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__weight_array__
#define __Feynman_Simulator__weight_array__

#include "utility/fft.h"
#include "utility/array.h"
#include "utility/convention.h"
#include "utility/complex.h"
#include "estimator/estimator.h"
#include "lattice/lattice.h"
#include <vector>

namespace weight0 {

const int MAX_BIN = 32;
typedef int Mode;
const Mode Spatial = 1;
const Mode Time = 2;
enum Dim { ORDER,
           SP,
           SUB,
           VOL,
           TAU
};

class WeightBase {
  protected:
    WeightBase(const Lattice &, real Beta, int SpinVol, const std::string &name);

    std::string _Name;
    real _Beta;
    real _dBeta; //_Beta/MAX_TAU
    real _dBetaInverse;
    Lattice _Lat;

    unsigned _Shape[5]; //the shape of internal weight array

    int SpinIndex(spin SpinIn, spin SpinOut);
    int SpinIndex(spin *TwoSpinIn, spin *TwoSpinOut);
    int TauToBin(real t);
    int TauToBin(real t_in, real t_out);
    real BinToTau(int Bin);
};

class WeightArray : public WeightBase {
  public:
    WeightArray(const Lattice &, real Beta, int SpinVol, const std::string &name);
    void SetTest();
    void Save(const std::string &FileName, std::string Mode = "a");
    bool Load(const std::string &);
    void Reset(real Beta);

    Array::array4<Complex> Weight;

    //No Time variable for this guy
    Array::array3<Complex> DeltaTWeight;

    //No Time variable for this guy as well
    //But you can not merge BareWeight with DeltaTWeight
    //You will have to do all kinds of manipulation(like fft) on _BareWeight
    Array::array3<Complex> BareWeight;

    unsigned int SpaceTimeShape[D + 1]; //store Lx,Ly,Lz,Lt

    unsigned int *Shape();
    unsigned int Size();

  protected:
    bool _CheckVec2Index();
    void _FFT(fft::Dir, Mode);
    void _ChangeSymmetry(fft::Dir);
};

class WeightEstimator : public WeightBase {
  public:
    WeightEstimator(const Lattice &, real Beta, int Spine, const std::string &name,
                    int order, real Norm);

    Estimate<Complex> WeightWithError(int order);

    int OrderAcceptable(int StartFromOrder, real ErrorThreshold);
    //return final weight density of tau
    void GetWeightArray(WeightArray &Weight, int UpToOrder);

    //The internal _Beta will be changed, so do _WeightAccu, _DeltaWeightAccu and _NormAccu
    //all changed will be done to make sure GetWeightArray returns the reweighted weight function
    //(as for now, reweighted weight function is set to be the unreweighted weight function)
    void ReWeight(real Beta);

    void AddStatistics();
    void ClearStatistics();
    void SqueezeStatistics(real factor);
    //    std::string PrettyString();
    void Save(const std::string &FileName, std::string Mode = "a");
    bool Load(const std::string &);

  protected:
    int _Order;
    //final weight of each bin = _WeightAccu/_NormAccu*_Norm
    //final weight function = (final weight of each bin)*MAX_BIN/Beta
    real _Norm;     //The normalization factor
    real _NormAccu; //The normalization accumulation
    Array::array5<Complex> _WeightAccu;

    unsigned int *Shape();
    EstimatorBundle<Complex> _Average;
};
}

#endif /* defined(__Feynman_Simulator__weight_array__) */
