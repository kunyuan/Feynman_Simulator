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
bool Convert(Object obj, std::string& val)
{
    if (!PyString_Check(obj.Borrow()))
        return false;
    val = PyString_AsString(obj.Borrow());
    return true;
}
bool Convert(Object obj, bool& value)
{
    if (obj.Borrow() == Py_False)
        value = false;
    else if (obj.Borrow() == Py_True)
        value = true;
    else
        return false;
    return true;
}
bool IsInt(const Object& obj)
{
    return (!PyInt_Check(obj.Borrow()) && !PyLong_Check(obj.Borrow()));
}

bool Convert(Object obj, int& value)
{
    if (IsInt(obj))
        return false;
    value = (int)PyLong_AsLong(obj.Borrow());
    return true;
}
bool Convert(Object obj, unsigned int& value)
{
    if (IsInt(obj))
        return false;
    value = (unsigned int)PyLong_AsUnsignedLong(obj.Borrow());
    return true;
}
bool Convert(Object obj, long& value)
{
    if (IsInt(obj))
        return false;
    value = PyLong_AsLong(obj.Borrow());
    return true;
}
bool Convert(Object obj, unsigned long& value)
{
    if (IsInt(obj))
        return false;
    value = PyLong_AsUnsignedLong(obj.Borrow());
    return true;
}
bool Convert(Object obj, long long& value)
{
    if (IsInt(obj))
        return false;
    value = PyLong_AsLongLong(obj.Borrow());
    return true;
}
bool Convert(Object obj, unsigned long long& value)
{
    if (IsInt(obj))
        return false;
    value = PyLong_AsUnsignedLongLong(obj.Borrow());
    return true;
}

template <typename T>
bool ConvertReal(Object obj, T& val)
{
    if (!PyFloat_Check(obj.Borrow()))
        return false;
    val = (T)PyFloat_AsDouble(obj.Borrow());
    return true;
}
bool Convert(Object obj, float& val)
{
    return ConvertReal(obj, val);
}
bool Convert(Object obj, double& val)
{
    return ConvertReal(obj, val);
}
bool Convert(Object obj, Complex& val)
{
    if (!PyComplex_Check(obj.Borrow()))
        return false;
    val = Complex(PyComplex_RealAsDouble(obj.Borrow()),
                  PyComplex_ImagAsDouble(obj.Borrow()));
    return true;
}

// Allocation methods

Object CastToPy(const std::string& str)
{
    return Object::Steal(PyString_FromString(str.c_str()));
}
Object CastToPy(bool value)
{
    return Object::Steal(PyBool_FromLong(value));
}
Object CastToPy(int num)
{
    return Object::Steal(PyInt_FromLong(num));
}
Object CastToPy(unsigned int num)
{
    return Object::Steal(PyInt_FromLong(num));
}
Object CastToPy(long num)
{
    return Object::Steal(PyLong_FromLong(num));
}
Object CastToPy(unsigned long num)
{
    return Object::Steal(PyLong_FromUnsignedLong(num));
}
Object CastToPy(unsigned long long num)
{
    return Object::Steal(PyLong_FromUnsignedLongLong(num));
}
Object CastToPy(long long num)
{
    return Object::Steal(PyLong_FromLongLong(num));
}
Object CastToPy(unsigned long long num);

Object CastToPy(float num)
{
    return Object::Steal(PyFloat_FromDouble(num));
}
Object CastToPy(double num)
{
    return Object::Steal(PyFloat_FromDouble(num));
}
Object CastToPy(const Complex& num)
{
    return Object::Steal(PyComplex_FromDoubles(num.Re, num.Im));
}
}
