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

class WeightNoMeasure {
    //dyson::Dyson does care abort the memroy structure of weight!!!
    friend dyson::Dyson;

  public:
    void SetTest();
    void Save(const std::string &FileName, std::string Mode = "a");
    bool Load(const std::string &);
    void Reset(real Beta);

    enum Dim { ORDER,
               SP,
               SUB,
               VOL,
               TAU };
    unsigned int *Shape();

  protected:
    WeightNoMeasure(const Lattice &, real Beta, int Order, int SpinVol, std::string);
    std::string _Name;
    real _Beta;
    real _dBeta; //_Beta/MAX_TAU
    real _dBetaInverse;
    int _Order;
    Lattice _Lat;
    Array::array4<Complex> _Weight;

    //No Time variable for this guy
    Array::array3<Complex> _DeltaTWeight;

    //No Time variable for this guy as well
    //But you can not merge _BareWeight with _DeltaTWeight
    //You will have to do all kinds of manipulation(like fft) on _BareWeight
    Array::array3<Complex> _BareWeight;

    unsigned int _Shape[5];              //the shape of internal weight array
    unsigned int _SpaceTimeShape[D + 1]; //store Lx,Ly,Lz,Lt
    bool _CheckVec2Index();
    void _FFT(fft::Dir, Mode);
    void _ChangeSymmetry(fft::Dir);

    int SpinIndex(spin SpinIn, spin SpinOut);
    int SpinIndex(spin *TwoSpinIn, spin *TwoSpinOut);
    int TauToBin(real tau);
    real BinToTau(int Bin);
};

class WeightNeedMeasure : public WeightNoMeasure {
  protected:
    real _Norm;
    Array::array5<Complex> _WeightAccu;
    EstimatorBundle<Complex> _Average;

  public:
    WeightNeedMeasure(const Lattice &, real Beta, int order, int Spine, std::string);

    Estimate<Complex> WeightWithError(int order);

    int OrderAcceptable(int StartFromOrder, real ErrorThreshold);
    void UpdateWeight(int UpToOrder);

    void MeasureNorm(); //weight=Beta/MAX_BIN/zeroth order weight
    real Normalization();

    void AddStatistics();
    void ClearStatistics();
    void SqueezeStatistics(real factor);
    //    std::string PrettyString();
    void Save(const std::string &FileName, std::string Mode = "a");
    bool Load(const std::string &);
};
}

#endif
