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
#include <math.h>

using namespace std;
using namespace weight;

const string SMOOTH = "SmoothT";
const string DELTA = "DeltaT";

bool DeltaTArray::FromDict(const Dictionary& dict)
{
    auto arr = dict.Get<Python::ArrayObject>(DELTA);
    ASSERT_ALLWAYS(Equal(arr.Shape().data(), GetShape(), GetDim()), "Shape should match!");
    Assign(arr.Data<Complex>());
    return true;
}
Dictionary DeltaTArray::ToDict()
{
    Dictionary dict;
    dict[DELTA] = Python::ArrayObject(Data(), GetShape(), GetDim());
    return dict;
}
bool SmoothTArray::FromDict(const Dictionary& dict)
{
    auto arr = dict.Get<Python::ArrayObject>(SMOOTH);
    ASSERT_ALLWAYS(Equal(arr.Shape().data(), GetShape(), GetDim()), "Shape should match!");
    Assign(arr.Data<Complex>());
    return true;
}
Dictionary SmoothTArray::ToDict()
{
    Dictionary dict;
    dict[SMOOTH] = Python::ArrayObject(Data(), GetShape(), GetDim());
    return dict;
}
