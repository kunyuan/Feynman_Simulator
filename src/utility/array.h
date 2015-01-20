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
#include <initializer_list>
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
    Array(const Array& source) = delete;
    Array& operator=(const Array& c) = delete;
    ~Array() { Free(); };
    void Allocate(const uint* Shape_);
    void Allocate(const std::vector<uint> Shape_)
    {
        Allocate(Shape_.data());
    }
    void Free();
    void Copy(const Array& c);
    void Assign(const Complex& c);
    void Assign(const Complex* c);
    void Assign(const std::initializer_list<Complex>& source);

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
            pos += index[i - 1] * _Cache[i];
        pos += index[DIM - 2];
        if (DEBUGMODE && pos >= _Size)
            THROW_ERROR(IndexInvalid, "exceed array bound!");
        return _Data[pos];
    }
    inline Complex& operator()(const uint* index)
    {
        uint pos = index[0] * _Cache[0];
        for (uint i = 1; i < DIM - 1; i++)
            pos += index[i] * _Cache[i];
        pos += index[DIM - 1];
        if (DEBUGMODE && pos >= _Size)
            THROW_ERROR(IndexInvalid, "exceed array bound!");
        return _Data[pos];
    }
    inline Complex At(const uint* index) const
    {
        uint pos = index[0] * _Cache[0];
        for (uint i = 1; i < DIM - 1; i++)
            pos += index[i] * _Cache[i];
        pos += index[DIM - 1];
        if (DEBUGMODE && pos >= _Size)
            THROW_ERROR(IndexInvalid, "exceed array bound!");
        return _Data[pos];
    }
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
int TestArray();
#endif /* defined(__Feynman_Simulator__array__) */
