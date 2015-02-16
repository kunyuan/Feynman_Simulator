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
#include <string>

class Dictionary;
namespace weight {

enum SpinNum {
    SPIN2 = 2,
    SPIN4 = 4
};

const std::string SMOOTH = "SmoothT";
const std::string DELTA = "DeltaT";

template <uint DIM>
class WeightArray {
public:
    WeightArray()
        : _Data(nullptr){};
    //copy sematics everywhere
    WeightArray(const WeightArray& source) = delete;
    WeightArray& operator=(const WeightArray& c) = delete;
    ~WeightArray() { Free(); };
    void Allocate(const uint* shape, const std::string Name);
    void Free();
    void Copy(const WeightArray& c);
    void Assign(const Complex& c);
    void Assign(const Complex* c); //copy _Size complex into _Data
    void Assign(const Complex* c, uint size); //copy size complex into _Data

    uint GetDim() const { return DIM; }
    uint GetSize() const { return _Size; }
    const uint* GetShape() const { return _Shape; }
    Complex& operator[](uint Index) { return _Data[Index]; }
    const Complex& operator()(uint Index) const { return _Data[Index]; }
    Complex* Data() { return _Data; }

    bool FromDict(const Dictionary&);
    Dictionary ToDict();

    template <typename T>
    WeightArray& operator+=(const T& rhs)
    {
        for (uint i = 0; i < _Size; i++)
            _Data[i] += rhs;
        return *this;
    }
    template <typename T>
    WeightArray& operator-=(const T& rhs)
    {
        for (uint i = 0; i < _Size; i++)
            _Data[i] -= rhs;
        return *this;
    }
    template <typename T>
    WeightArray& operator*=(const T& rhs)
    {
        for (uint i = 0; i < _Size; i++)
            _Data[i] *= rhs;
        return *this;
    }
    template <typename T>
    WeightArray& operator/=(const T& rhs)
    {
        for (uint i = 0; i < _Size; i++)
            _Data[i] /= rhs;
        return *this;
    }

protected:
    Complex* _Data;
    bool IsAllocated;
    uint _Shape[DIM];
    uint _Size;
    std::string _Name;
};
}

#endif /* defined(__Feynman_Simulator__weight_basic__) */
