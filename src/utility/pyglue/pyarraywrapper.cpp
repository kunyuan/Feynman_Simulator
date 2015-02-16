//
//  pyarraywrapper.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 12/25/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "utility/complex.h"
#include "pyarraywrapper.h"
#include <Python.h>
#include <numpy/arrayobject.h>

using namespace std;
namespace Python {

void ArrayInitialize()
{
    if (_import_array() < 0) {
        PropagatePyError();
        PyErr_Print();
        PyErr_SetString(PyExc_ImportError, "numpy.core.multiarray failed to import");
    }
}

bool Convert(Object obj, ArrayObject& array)
{
    //    if (array_init == 0)
    if (!PyArray_Check(obj.Get()))
        return false;
    array = obj;
    return true;
}

ArrayObject::ArrayObject(const Object& obj)
{
    if (!PyArray_Check(obj.Get()))
        ABORT("PyArray object is expected!");
    if (PyArray_IS_C_CONTIGUOUS(obj.Get()))
        _PyPtr = obj.Get(NewRef);
    else
        //NewCopy return a new reference
        _PyPtr = PyArray_NewCopy((PyArrayObject*)obj.Get(), NPY_CORDER);
}
ArrayObject::ArrayObject(PyObject* obj, OwnerShip ownership)
{
    Object temp(obj, ownership);
    if (!PyArray_Check(obj))
        ABORT("PyArray object is expected!");
    if (PyArray_IS_C_CONTIGUOUS(obj))
        _PyPtr = temp.Get(NewRef);
    else
        //NewCopy return a new reference
        _PyPtr = PyArray_NewCopy((PyArrayObject*)temp.Get(), NPY_CORDER);
}

void ArrayObject::_Construct(Complex* data, const uint* Shape, const int Dim)
{
    ASSERT_ALLWAYS(data != nullptr, "data pointer shouldn't be null!");
    ASSERT_ALLWAYS(Shape != nullptr, "Shape pointer shouldn't be null!");
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
    //A new reference to an ndarray is returned, but the ndarray will not own its data.
    //When this ndarray is deallocated, the pointer will not be freed.
    //You should ensure that the provided memory is not freed while the returned array is in existence.
    PropagatePyError();
    ASSERT_ALLWAYS(array != nullptr, "Failed to create python array!");
    *this = Object(array);
}

void ArrayObject::_Construct(real* data, const uint* Shape, const int Dim)
{
    ASSERT_ALLWAYS(data != nullptr, "data pointer shouldn't be null!");
    ASSERT_ALLWAYS(Shape != nullptr, "Shape pointer shouldn't be null!");
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
    //A new reference to an ndarray is returned, but the ndarray will not own its data.
    //When this ndarray is deallocated, the pointer will not be freed.
    //You should ensure that the provided memory is not freed while the returned array is in existence.
    PropagatePyError();
    ASSERT_ALLWAYS(array != nullptr, "Failed to create python array!");
    *this = Object(array);
}

template <>
Complex* ArrayObject::Data<Complex>()
{
    ASSERT_ALLWAYS(_PyPtr != nullptr, "ArrayObject is still empty!");
    return reinterpret_cast<Complex*>(PyArray_DATA(_PyPtr));
}
template <>
real* ArrayObject::Data<real>()
{
    ASSERT_ALLWAYS(_PyPtr != nullptr, "ArrayObject is still empty!");
    return reinterpret_cast<real*>(PyArray_DATA(_PyPtr));
}

std::vector<uint> ArrayObject::Shape()
{
    vector<uint> _Shape;
    ASSERT_ALLWAYS(_PyPtr != nullptr, "ArrayObject is still empty!");
    auto dim = PyArray_DIMS(_PyPtr);
    for (int i = 0; i < Dim(); i++)
        _Shape.push_back((uint)dim[i]);
    return _Shape;
}

uint ArrayObject::Size()
{
    uint size = 1;
    for (auto i : Shape())
        size *= i;
    return size;
}

int ArrayObject::Dim()
{
    ASSERT_ALLWAYS(_PyPtr != nullptr, "ArrayObject is still empty!");
    return PyArray_NDIM(_PyPtr);
}
}
