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
    //copy sematics everywhere
    Array(const Array& source);
    Array& operator=(const Complex& c);
    Array& operator=(const Array& c);
    Array& operator=(const Complex* c);
    ~Array() { Free(); };

    const uint* GetShape() const;
    Complex* Data() const;
    Complex& operator()(const std::vector<uint> index)
    {
        return (*this)(index.data());
    }
    inline Complex& operator()(uint index0, const uint* index)
    {
        uint pos = index0 * _Cache[0];
        for (uint i = 1; i < DIM - 1; i++)
            pos += index[i] * _Cache[i];
        return _Data[pos + index[DIM - 1]];
    }
    inline Complex& operator()(const uint* index)
    {
        uint pos = index[0] * _Cache[0];
        for (uint i = 1; i < DIM - 1; i++)
            pos += index[i] * _Cache[i];
        return _Data[pos + index[DIM - 1]];
    }
    inline Complex At(const uint* index) const
    {
        uint pos = index[0] * _Cache[0];
        for (uint i = 1; i < DIM - 1; i++)
            pos += index[i] * _Cache[i];
        return _Data[pos + index[DIM - 1]];
    }
    void Allocate(const uint* Shape_, const Complex* data = nullptr);
    void Free();

    template <typename T>
    Array& operator+=(const T& rhs)
    {
        for (uint i = 0; i < _Size; i++)
            _Data[i] += rhs;
        return *this;
    }
    template <typename T>
    Array& operator-=(const T& rhs)
    {
        for (uint i = 0; i < _Size; i++)
            _Data[i] -= rhs;
        return *this;
    }
    template <typename T>
    Array& operator*=(const T& rhs)
    {
        for (uint i = 0; i < _Size; i++)
            _Data[i] *= rhs;
        return *this;
    }
    template <typename T>
    Array& operator/=(const T& rhs)
    {
        for (uint i = 0; i < _Size; i++)
            _Data[i] /= rhs;
        return *this;
    }
};
#endif /* defined(__Feynman_Simulator__array__) */
