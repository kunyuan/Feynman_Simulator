//
//  pyarraywrapper.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 12/25/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__pyarraywrapper__
#define __Feynman_Simulator__pyarraywrapper__

#include "utility/abort.h"
#include <vector>
#include "object.h"

namespace Python {
void ArrayInitialize();
class ArrayObject;
bool Convert(Object, ArrayObject&);
//Object CastToPy(const ArrayObject&);
class ArrayObject : public Object {
public:
    ArrayObject()
        : Object()
    {
    }
    ArrayObject(const Object& obj);
    ArrayObject(PyObject* obj, OwnerShip ownership = NewRef);
    template <typename T>
    ArrayObject(T* data, const std::vector<uint>& Shape, const int Dim)
    {
        _Construct(data, Shape.data(), Dim);
    }
    template <typename T>
    ArrayObject(T* data, const uint* Shape, const int Dim)
    {
        _Construct(data, Shape, Dim);
    }
    template <typename T>
    T* Data();
    std::vector<uint> Shape();
    uint Size();
    int Dim();
    ArrayObject& operator=(const ArrayObject& obj)
    {
        Object::operator=(obj);
        return *this;
    }

private:
    void _Construct(real* data, const uint* Shape, const int Dim);
    void _Construct(Complex* data, const uint* Shape, const int Dim);
};
}

#endif /* defined(__Feynman_Simulator__pyarraywrapper__) */
