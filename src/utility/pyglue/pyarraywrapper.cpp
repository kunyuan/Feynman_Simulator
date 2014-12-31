//
//  pyarraywrapper.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 12/25/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "utility/complex.h"
#include "pyarraywrapper.h"
#include <Python/Python.h>
#include <numpy/arrayobject.h>

using namespace std;
namespace Python {

void ArrayInitialize()
{
    import_array();
}

bool Convert(Object obj, ArrayObject& array)
{
    if (!PyArray_Check(obj.Get()))
        return false;
    array = obj;
    return true;
}

ArrayObject::ArrayObject(const Object& obj)
    : Object(obj)
{
    if (!PyArray_Check(obj.Get()))
        ERRORCODEABORT(ERR_VALUE_INVALID, "PyArray object is expected!");
}
ArrayObject::ArrayObject(PyObject* obj, OwnerShip ownership)
    : Object(obj, ownership)
{
    if (!PyArray_Check(obj))
        ERRORCODEABORT(ERR_VALUE_INVALID, "PyArray object is expected!");
}

void ArrayObject::_Construct(Complex* data, const uint* Shape, const int Dim)
{
    int TypeName;
    if (sizeof(data[0]) == 8)
        TypeName = NPY_COMPLEX64;
    else if (sizeof(data[0]) == 16)
        TypeName = NPY_COMPLEX128;
    else
        ABORT("why complex number has size>16?");
    vector<npy_intp> _Shape;
    for (int i = 0; i < Dim; i++)
        _Shape.push_back((npy_intp)Shape[i]);
    PyObject* array = PyArray_SimpleNewFromData(Dim, _Shape.data(), TypeName, (void*)data);
    *this = Object(array);
}

void ArrayObject::_Construct(real* data, const uint* Shape, const int Dim)
{
    int TypeName;
    if (sizeof(data[0]) == 4)
        TypeName = NPY_FLOAT32;
    else if (sizeof(data[0]) == 8)
        TypeName = NPY_FLOAT64;
    else
        ABORT("why real number has size>8?");
    vector<npy_intp> _Shape;
    for (int i = 0; i < Dim; i++)
        _Shape.push_back((npy_intp)Shape[i]);
    PyObject* array = PyArray_SimpleNewFromData(Dim, _Shape.data(), TypeName, (void*)data);
    *this = Object(array);
}

template <>
Complex* ArrayObject::Data<Complex>()
{
    return reinterpret_cast<Complex*>(PyArray_DATA(_PyPtr));
}
template <>
real* ArrayObject::Data<real>()
{
    return reinterpret_cast<real*>(PyArray_DATA(_PyPtr));
}

std::vector<uint> ArrayObject::Shape()
{
    vector<uint> _Shape;
    auto dim = PyArray_DIMS(_PyPtr);
    for (int i = 0; i < Dim(); i++)
        _Shape.push_back((uint)dim[i]);
    return _Shape;
}

int ArrayObject::Dim()
{
    return PyArray_NDIM(_PyPtr);
}
}
