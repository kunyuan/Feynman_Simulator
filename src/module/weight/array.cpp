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
uint* Array<DIM>::GetShape()
{
    return _Shape;
}

template <uint DIM>
Complex* Array<DIM>::Data()
{
    return _Data;
}

template <uint DIM>
Complex& Array<DIM>::operator()(uint* index)
{
    uint pos;
    for (uint i = 0; i < DIM - 1; i++)
        pos += index[i] * _Cache[i];
    return _Data[pos + index[DIM - 1]];
}

template <uint DIM>
void Array<DIM>::Allocate(uint* Shape_, Complex* data)
{
    ASSERT_ALLWAYS(_Data == nullptr, "Please free Array first!");
    std::copy(Shape_, Shape_ + DIM, _Shape);
    _Size = 1;
    for (uint i = 0; i < DIM; i++) {
        _Cache[DIM - 1 - i] = _Size;
        _Size *= _Shape[i];
    }
    _Data = new Complex[_Size];
    if(_Data==nullptr)
        THROW_ERROR(MemoryException, "Fail to allocate array!");
    std::copy(data, data + DIM, _Data);
}

template <uint DIM>
void Array<DIM>::Free()
{
    delete[] _Data;
    _Data = nullptr;
}
