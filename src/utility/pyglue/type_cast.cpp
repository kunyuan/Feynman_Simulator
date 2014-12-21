//
//  type_cast.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 12/20/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "type_cast.h"
#include "utility/complex.h"
namespace Python {
bool Convert(PyObject* obj, std::string& val)
{
    if (!PyString_Check(obj))
        return false;
    val = PyString_AsString(obj);
    return true;
}
bool Convert(PyObject* obj, bool& value)
{
    if (obj == Py_False)
        value = false;
    else if (obj == Py_True)
        value = true;
    else
        return false;
    return true;
}
bool Convert(PyObject* obj, int& value)
{
    if (!PyInt_Check(obj) && !PyLong_Check(obj))
        return false;
    value = (int)PyLong_AsLong(obj);
    return true;
}
bool Convert(PyObject* obj, unsigned int& value)
{
    if (!PyInt_Check(obj) && !PyLong_Check(obj))
        return false;
    value = (unsigned int)PyLong_AsUnsignedLong(obj);
    return true;
}
bool Convert(PyObject* obj, long& value)
{
    if (!PyInt_Check(obj) && !PyLong_Check(obj))
        return false;
    value = PyLong_AsLong(obj);
    return true;
}
bool Convert(PyObject* obj, unsigned long& value)
{
    if (!PyInt_Check(obj) && !PyLong_Check(obj))
        return false;
    value = PyLong_AsUnsignedLong(obj);
    return true;
}
bool Convert(PyObject* obj, long long& value)
{
    if (!PyInt_Check(obj) && !PyLong_Check(obj))
        return false;
    value = PyLong_AsLongLong(obj);
    return true;
}
bool Convert(PyObject* obj, unsigned long long& value)
{
    if (!PyInt_Check(obj) && !PyLong_Check(obj))
        return false;
    value = PyLong_AsUnsignedLongLong(obj);
    return true;
}

template <typename T>
bool ConvertReal(PyObject* obj, T& val)
{
    if (!PyFloat_Check(obj))
        return false;
    val = (T)PyFloat_AsDouble(obj);
    return true;
}
bool Convert(PyObject* obj, float& val)
{
    return ConvertReal(obj, val);
}
bool Convert(PyObject* obj, double& val)
{
    return ConvertReal(obj, val);
}
bool Convert(PyObject* obj, Complex& val)
{
    if (!PyComplex_Check(obj))
        return false;
    val = Complex(PyComplex_RealAsDouble(obj),
                  PyComplex_ImagAsDouble(obj));
    return true;
}

// Allocation methods

PyObject* CastToPyObject(const std::string& str)
{
    return PyString_FromString(str.c_str());
}
PyObject* CastToPyObject(bool value)
{
    return PyBool_FromLong(value);
}
PyObject* CastToPyObject(int num)
{
    return PyInt_FromLong(num);
}
PyObject* CastToPyObject(unsigned int num)
{
    return PyInt_FromLong(num);
}
PyObject* CastToPyObject(long num)
{
    return PyLong_FromLong(num);
}
PyObject* CastToPyObject(unsigned long num)
{
    return PyLong_FromUnsignedLong(num);
}
PyObject* CastToPyObject(unsigned long long num)
{
    return PyLong_FromUnsignedLongLong(num);
}
PyObject* CastToPyObject(long long num)
{
    return PyLong_FromLongLong(num);
}
PyObject* CastToPyObject(unsigned long long num);

PyObject* CastToPyObject(float num)
{
    return PyFloat_FromDouble(num);
}
PyObject* CastToPyObject(double num)
{
    return PyFloat_FromDouble(num);
}
PyObject* CastToPyObject(const Complex& num)
{
    return PyComplex_FromDoubles(num.Re, num.Im);
}
}
