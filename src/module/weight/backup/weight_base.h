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
#include <vector>

namespace dyson {
class Dyson;
}

namespace weight {

const int MAX_BIN = 32;
typedef int Mode;
const Mode Spatial = 1;
const Mode Time = 2;

enum Dim {
    SP,
    SUB,
    VOL,
    TAU,
};

class WeightNoMeasure {
  public:
    void Save(const std::string &FileName, std::string Mode = "a");
    bool Load(const std::string &);
    void Reset(real Beta);
    unsigned int *Shape();

    Array::array4<Complex> SmoothWeight;
    //No Time variable for this guy
    Array::array3<Complex> DeltaTWeight;

    int SpinIndex(spin SpinIn, spin SpinOut);
    bool IsSameSpin(int spindex)
    {
        return (spindex == 0 || spindex == 2);
    }

    //First In/Out: direction of WLine; Second In/Out: direction of Vertex
    int SpinIndex(spin SpinInIn, spin SpinInOut, spin SpinOutIn, spin SpinOutOut);
    int SpinIndex(spin *TwoSpinIn, spin *TwoSpinOut);

    enum SpinFilter { UpUp2UpUp,
                      UpDown2UpDown,
                      UpDown2DownUp };
    std::vector<int> GetSpinIndexVector_FourSpinsFileter(SpinFilter);

  protected:
    WeightNoMeasure(const Lattice &, real Beta,
                    bool IsTauSymmetric, int SpinVol, std::string);
    model _Model;
    bool _IsTauSymmetric;
    std::string _Name;
    real _Beta;
    real _dBeta; //_Beta/MAX_TAU
    real _dBetaInverse;
    Lattice _Lat;

    unsigned int _Shape[4];              //the shape of internal weight array
    unsigned int _SpaceTimeShape[D + 1]; //store Lx,Ly,Lz,Lt
    bool _CheckVec2Index();
    void _FFT(fft::Dir, Mode);
    void _ChangeSymmetry(fft::Dir);

    /**
    *  @return return 1 if weight is symmtric in tau, otherwise, return -1
    */
    int TauSymmetry(real t_in, real t_out);
    int TauToBin(real tau);
    int TauToBin(real t_in, real t_out);
    real BinToTau(int Bin);
};
}

#endif
