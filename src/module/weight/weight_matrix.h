//
//  weight_array.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/20/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__weight_array__
#define __Feynman_Simulator__weight_array__

#include "utility/array.h"
#include "utility/complex.h"
#include "utility/fft.h"
#include "weight_basic.h"
#include <vector>

class Lattice;
namespace weight {
const uint MAX_TAU_BIN = 32;
enum FFT_Mode {
    Spatial = 1,
    Time = 2
};
enum SpinNum {
    OneSpin = 1,
    TwoSpins = 2
};
enum Dim {
    SP,
    SUB,
    VOL,
    TAU,
};

/**
*  This is the weight function, not array of numbers, so you should take care of dBeta inside
*/
class SmoothTMatrix : public Array::array4<Complex> {
  public:
    /**
    *  SoomthTWeight will not allocate memory!!!
    */
    SmoothTMatrix(const Lattice &lat, SpinNum, const std::string Name);
    /**
    *  SoomthTWeight will allocate memory after calling Activate(), but not initialize elements!!!
    */
    void Activate();
    vector<uint> Shape;
    vector<uint> SpaceTimeShape;
    void FFT(FFT_Mode, fft::Dir);
    void MatrixInverse();

    bool Load(const std::string &FileName);
    void Save(const std::string &FileName, const std::string Mode = "a");

  private:
    std::string _Name;
};

class DeltaTMatrix : public Array::array3<Complex> {
  public:
    DeltaTMatrix(const Lattice &lat, SpinNum, const std::string Name);
    void Activate();
    vector<uint> Shape;
    vector<uint> SpaceShape;
    void FFT(FFT_Mode, fft::Dir);

    bool Load(const std::string &FileName);
    void Save(const std::string &FileName, const std::string Mode = "a");

  private:
    std::string _Name;
};
}

#endif /* defined(__Feynman_Simulator__weight_array__) */
