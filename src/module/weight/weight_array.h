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

enum SpinNum {
    SPIN2 = 2,
    SPIN4 = 4
};

enum Dim {
    SP1 = 0,
    SUB1,
    SP2,
    SUB2,
    VOL,
    TAU,
};

class DeltaTArray : public Array<DELTA_T_SIZE> {
public:
    DeltaTArray()
        : Array<DELTA_T_SIZE>()
    {
    }
    DeltaTArray(const DeltaTArray&) = delete;
    DeltaTArray& operator=(const DeltaTArray& c) = delete;
    bool FromDict(const Dictionary&);
    Dictionary ToDict();
};
class SmoothTArray : public Array<SMOOTH_T_SIZE> {
public:
    SmoothTArray()
        : Array<SMOOTH_T_SIZE>()
    {
    }
    SmoothTArray(const SmoothTArray&) = delete;
    SmoothTArray& operator=(const SmoothTArray& c) = delete;
    bool FromDict(const Dictionary&);
    Dictionary ToDict();
};
}

#endif /* defined(__Feynman_Simulator__weight_basic__) */
