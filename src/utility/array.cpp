//
//  array.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 1/18/15.
//  Copyright (c) 2015 Kun Chen. All rights reserved.
//

#include <algorithm>
#include "array.h"
#include "utility/abort.h"
#include "utility/dictionary.h"

template <uint DIM>
void Array<DIM>::Copy(const Array& source)
{
    if (&source == this)
        return;
    if (source.Data() == nullptr) {
        Free();
        return;
    }
    if (_Data == nullptr)
        Allocate(source.GetShape());
    else
        ASSERT_ALLWAYS(Equal(_Shape, source.GetShape(), DIM), "Shape should be allocated first!");
    Assign(source.Data());
}
template <uint DIM>
void Array<DIM>::Assign(const Complex& c)
{
    ASSERT_ALLWAYS(_Data != nullptr, "Array should be allocated first!");
    for (uint i = 0; i < _Size; i++)
        _Data[i] = c;
}
template <uint DIM>
void Array<DIM>::Assign(const Complex* source)
{
    ASSERT_ALLWAYS(_Data != nullptr, "Array should be allocated first!");
    if (_Data == source)
        return;
    std::copy(source, source + _Size, _Data);
}

template <uint DIM>
const uint* Array<DIM>::GetShape() const
{
    return _Shape;
}

template <uint DIM>
Complex* Array<DIM>::Data() const
{
    return _Data;
}

template <uint DIM>
void Array<DIM>::Allocate(const uint* Shape_)
{
    ASSERT_ALLWAYS(_Data == nullptr, "Please free Array first!");
    std::copy(Shape_, Shape_ + DIM, _Shape);
    _Size = 1;
    for (auto i = 0; i < DIM; i++) {
        _Cache[DIM - 1 - i] = _Size;
        _Size *= _Shape[DIM - 1 - i];
    }
    _Data = new Complex[_Size];
    if (_Data == nullptr)
        THROW_ERROR(MemoryException, "Fail to allocate array!");
}

template <uint DIM>
void Array<DIM>::Free()
{
    delete[] _Data;
    _Data = nullptr;
}

template class Array<1>;
template class Array<2>;
template class Array<3>;
template class Array<4>;
template class Array<5>;
