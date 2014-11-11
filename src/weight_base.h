//
//  weight_base.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/7/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef Feynman_Simulator_weight_base_h
#define Feynman_Simulator_weight_base_h
#include "utility/convention.h"
#include "utility/array.h"
#include "lattice/lattice.h"
#include "utility/complex.h"
#include "utility/fft.h"

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
  public:
    void SetTest();
    void Save(const std::string &FileName, std::string Mode = "a");
    bool Load(const std::string &);
    void Reset(real Beta);
    unsigned int *Shape();

    Array::array4<Complex> SmoothWeight;

    //No Time variable for this guy
    Array::array3<Complex> DeltaTWeight;


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
    int TauSymmetry(real t_in, real t_out);
    int TauToBin(real tau);
    int TauToBin(real t_in, real t_out);
    real BinToTau(int Bin);
};
}

#endif
