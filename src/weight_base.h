//
//  weight_base.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/7/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef Feynman_Simulator_weight_base_h
#define Feynman_Simulator_weight_base_h
namespace dyson {
class Dyson;
}

namespace weight {

const int MAX_BIN = 32;
typedef int Mode;
const Mode Spatial = 1;
const Mode Time = 2;

enum Dim { ORDER,
           SP,
           SUB,
           VOL,
           TAU };

class WeightNoMeasure {
    //dyson::Dyson does care abort the memroy structure of weight!!!
    friend dyson::Dyson;

  public:
    void SetTest();
    void Save(const std::string &FileName, std::string Mode = "a");
    bool Load(const std::string &);
    void Reset(real Beta);
    unsigned int *Shape();

    Array::array4<Complex> SmoothWeight;

    //No Time variable for this guy
    Array::array3<Complex> DeltaTWeight;

    //No Time variable for this guy as well
    //But you can not merge _BareWeight with _DeltaTWeight
    //You will have to do all kinds of manipulation(like fft) on _BareWeight
    Array::array3<Complex> BareWeight;

  protected:
    WeightNoMeasure(const Lattice &, real Beta, int Order, int SpinVol, std::string);
    std::string _Name;
    real _Beta;
    real _dBeta; //_Beta/MAX_TAU
    real _dBetaInverse;
    int _Order;
    Lattice _Lat;

    unsigned int _Shape[5];              //the shape of internal weight array
    unsigned int _SpaceTimeShape[D + 1]; //store Lx,Ly,Lz,Lt
    bool _CheckVec2Index();
    void _FFT(fft::Dir, Mode);
    void _ChangeSymmetry(fft::Dir);

    int SpinIndex(spin SpinIn, spin SpinOut);
    int SpinIndex(spin *TwoSpinIn, spin *TwoSpinOut);
    int TauToBin(real tau);
    int TauToBin(real t_in, real t_out);
    real BinToTau(int Bin);
};

class WeightNeedMeasure : public WeightNoMeasure {
  public:
    WeightNeedMeasure(const Lattice &, real Beta, int order,
                      int Spine, std::string, real Norm);

    unsigned int *Shape();

    Estimate<Complex> WeightWithError(int order);

    int OrderAcceptable(int StartFromOrder, real ErrorThreshold);
    //update final weight density to WeightNoMeasure._Weight
    void UpdateWeight(int UpToOrder);

    //The internal _Beta will be changed, so do _WeightAccu, _DeltaWeightAccu and _NormAccu
    //all changed will be done to make sure GetWeightArray returns the reweighted weight function
    //(as for now, reweighted weight function is set to be the unreweighted weight function)
    void ReWeight(real Beta);

    void MeasureNorm();
    void AddStatistics();
    void ClearStatistics();
    void SqueezeStatistics(real factor);
    //    std::string PrettyString();
    void Save(const std::string &FileName, std::string Mode = "a");
    bool Load(const std::string &);

  protected:
    //final weight of each bin = _WeightAccu/_NormAccu*_Norm
    //final weight function = (final weight of each bin)*MAX_BIN/Beta
    real _Norm;     //The normalization factor
    real _NormAccu; //The normalization accumulation
    Array::array5<Complex> _WeightAccu;
    EstimatorBundle<Complex> _Average;
};
}

#endif
