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

class Dictionary;
namespace weight {

enum TauSymmetry {
    TauSymmetric = 1,
    TauAntiSymmetric = -1
};

enum SpinNum {
    SPIN2 = 2,
    SPIN4 = 4
};

enum Dim {
    SP,
    SUB,
    VOL,
    TAU,
};
typedef Array<3> DeltaTMatrix;
typedef Array<4> SmoothTMatrix;
class Basic {
protected:
    Basic(const Lattice& lat, real Beta, uint MaxTauBin, SpinNum,
          TauSymmetry Symmetry, std::string);

    uint* GetShape(); //the shape of internal weight array
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
};
class DeltaTArray : public Array<3> {
public:
    DeltaTArray()
        : Array<3>()
    {
    }
    DeltaTArray(const DeltaTArray&) = delete;
    DeltaTArray& operator=(const DeltaTArray& c) = delete;
    bool FromDict(const Dictionary&);
    Dictionary ToDict();
};
class SmoothTArray : public Array<4> {
public:
    SmoothTArray()
        : Array<4>()
    {
    }
    SmoothTArray(const SmoothTArray&) = delete;
    SmoothTArray& operator=(const SmoothTArray& c) = delete;
    bool FromDict(const Dictionary&);
    Dictionary ToDict();
};
}

#endif /* defined(__Feynman_Simulator__weight_basic__) */
