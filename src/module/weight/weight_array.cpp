//
//  weight_basic.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/21/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "weight_array.h"
#include "utility/logger.h"
#include "utility/dictionary.h"
#include "index_map.h"
#include <math.h>

using namespace std;

namespace weight {
template <uint DIM>
void WeightArray<DIM>::Assign(const Complex& c)
{
    ASSERT_ALLWAYS(IsAllocated, "Array should be allocated first!");
    for (uint i = 0; i < _Size; i++)
        _Data[i] = c;
}
template <uint DIM>
void WeightArray<DIM>::Assign(const Complex* source)
{
    ASSERT_ALLWAYS(IsAllocated, "Array should be allocated first!");
    if (_Data == source)
        return;
    std::copy(source, source + _Size, _Data);
}

template <uint DIM>
void WeightArray<DIM>::Assign(const Complex* source, uint size)
{
    ASSERT_ALLWAYS(IsAllocated, "Array should be allocated first!");
    if (_Data == source)
        return;
    std::copy(source, source + size, _Data);
}

template <uint DIM>
void WeightArray<DIM>::Allocate(const uint* Shape_, const std::string Name)
{
    _Name = Name;
    if (IsAllocated)
        Free();
    std::copy(Shape_, Shape_ + DIM, _Shape);
    _Size = 1;
    for (auto i = 0; i < DIM; i++) {
        _Size *= _Shape[i];
    }
    _Data = new Complex[_Size];
    if (_Data == nullptr) {
        THROW_ERROR(MemoryException, "Fail to allocate array!");
        IsAllocated = false;
    }
    IsAllocated = true;
}

template <uint DIM>
void WeightArray<DIM>::Free()
{
    if (IsAllocated) {
        delete[] _Data;
        IsAllocated = false;
    }
}

template <uint DIM>
bool WeightArray<DIM>::FromDict(const Dictionary& dict)
{
    ASSERT_ALLWAYS(IsAllocated, "Array should be allocated first!");
    Python::ArrayObject arr = dict.Get<Python::ArrayObject>(_Name);
    ASSERT_ALLWAYS(Equal(arr.Shape().data(), GetShape(), GetDim()), "Shape should match!");
    Assign(arr.Data<Complex>());
    return true;
}

template <uint DIM>
Dictionary WeightArray<DIM>::ToDict()
{
    Dictionary dict;
    dict[_Name] = Python::ArrayObject(_Data, GetShape(), GetDim());
    return dict;
}

template class WeightArray<DELTA_T_SIZE>;
template class WeightArray<SMOOTH_T_SIZE>;
template class WeightArray<SMOOTH_T_SIZE + 1>;
}