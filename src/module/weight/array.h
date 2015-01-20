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
class Dictionary;
template <uint DIM>
class Array {
private:
    uint _Cache[DIM];
    uint _Shape[DIM];
    Complex* _Data;
    uint _Size;

public:
    Array()
        : _Data(nullptr){};
    ~Array() { Free(); };
    uint* GetShape();
    Complex* Data();
    Complex& operator()(uint* index);
    void Allocate(uint* Shape_, Complex* data = nullptr);
    void Free();
};
#endif /* defined(__Feynman_Simulator__array__) */
