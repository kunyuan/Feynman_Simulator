//
//  weight_basic.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/21/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__weight_basic__
#define __Feynman_Simulator__weight_basic__

#include "utility/convention.h"
#include "lattice/lattice.h"
#include "utility/complex.h"
#include "utility/fft.h"
#include <vector>
namespace weight0 {

const uint MAX_TAU_BIN = 32;

enum TauSymmetry {
    TauSymmetric = 1,
    TauAntiSymmetric = -1
};

enum SpinNum {
    TwoSpins = 2,
    FourSpins = 4
};

enum FFT_Mode {
    Spatial = 1,
    Time = 2
};

enum Dim {
    SP,
    SUB,
    VOL,
    TAU,
};

class Basic {
  protected:
    Basic(const Lattice &lat, real Beta, SpinNum, model Model,
          TauSymmetry Symmetry, std::string);

    vector<uint> GetShape();          //the shape of internal weight array
    vector<uint> GetSpaceShape();     //store Lx,Ly,Lz
    vector<uint> GetSpaceTimeShape(); //store Lx,Ly,Lz,Lt
    void Reset(real beta);

    int TauIndex(real tau);
    int TauIndex(real t_in, real t_out);
    real IndexToTau(int TauIndex);

    int SublatIndex(const Distance &dist);
    int CoordiIndex(const Distance &dist);

    model _Model;
    int _TauSymmetryFactor;
    std::string _Name;
    real _Beta;
    real _dBeta; //_Beta/MAX_TAU
    real _dBetaInverse;
    Lattice _Lat;
    int _SpinNum;
};

class BasicWithTwoSpins : public Basic {
  protected:
    BasicWithTwoSpins(const Lattice &lat, real Beta, model Model,
                      TauSymmetry Symmetry, std::string Name);
    static int SpinIndex(spin SpinIn, spin SpinOut);
    static bool IsSameSpin(int spindex);
};

class BasicWithFourSpins : public Basic {
  protected:
    BasicWithFourSpins(const Lattice &lat, real Beta, model Model,
                       TauSymmetry Symmetry, std::string Name);
    //First In/Out: direction of WLine; Second In/Out: direction of Vertex
    static int SpinIndex(spin SpinInIn, spin SpinInOut, spin SpinOutIn, spin SpinOutOut);
    static int SpinIndex(spin *TwoSpinIn, spin *TwoSpinOut);
    enum SpinFilter { UpUp2UpUp,
                      UpDown2UpDown,
                      UpDown2DownUp };
    static std::vector<int> GetSpinIndexVector(SpinFilter filter);
};
}

#endif /* defined(__Feynman_Simulator__weight_basic__) */
