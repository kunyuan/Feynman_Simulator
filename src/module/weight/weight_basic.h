//
//  weight_basic.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/21/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__weight_basic__
#define __Feynman_Simulator__weight_basic__

#include "utility/complex.h"
#include "utility/array.h"
#include "index_map.h"

namespace weight {

enum TauSymmetry {
    TauSymmetric = 1,
    TauAntiSymmetric = -1
};

enum SpinNum {
    SPIN2 = 2,
    SPIN4 = 4
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
typedef Array::array3<Complex> DeltaTMatrix;
typedef Array::array4<Complex> SmoothTMatrix;

class Basic {
    friend class GInitializer;
    friend class WInitializer;

public:
    bool Load(const std::string& FileName);
    void Save(const std::string& FileName, const std::string Mode = "a");

protected:
    Basic(const Lattice& lat, real Beta, uint MaxTauBin, SpinNum,
          TauSymmetry Symmetry, std::string);

    uint* GetShape(); //the shape of internal weight array
    uint* GetSpaceShape(); //store Lx,Ly,Lz
    uint* GetSpaceTimeShape(); //store Lx,Ly,Lz,Lt
    void Reset(real beta);
    int GetTauSymmetryFactor(real t_in, real t_out) const;

    int _TauSymmetryFactor;
    std::string _Name;
    real _Beta;
    real _dBeta; //_Beta/MAX_TAU
    real _dBetaInverse;
    uint _MaxTauBin;
    Lattice _Lat;
    int _SpinNum;
    vector<uint> _Shape;
    vector<uint> _SpaceTimeShape;

    DeltaTMatrix _DeltaTWeight;
    SmoothTMatrix _SmoothTWeight;
};
}

#endif /* defined(__Feynman_Simulator__weight_basic__) */
