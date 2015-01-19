//
//  array.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 1/18/15.
//  Copyright (c) 2015 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__array__
#define __Feynman_Simulator__array__
#include "utility/complex.h"
#include <vector>
class Dictionary;
namespace weight {
enum SpinNum {
    SPIN2 = 2,
    SPIN4 = 4
};

template <uint D>
class Array {
public:
    Complex* Data;
    uint Size;
    uint Dim;
    uint* Shape;
    void BuildNew(const std::vector<uint>& Shape);
    bool FromDict(const Dictionary&);
    Dictionary ToDict();
};
//
//class DeltaTArray : public Array {
//public:
//    DeltaTArray(SpinNum, IndexMap);
//    bool FromDict(const Dictionary&);
//    Dictionary ToDict();
//};
//
//class SmoothTArray : public Array {
//};
}
#endif /* defined(__Feynman_Simulator__array__) */
