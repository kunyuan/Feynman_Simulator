//
//  type_cast.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 12/20/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "utility/complex.h"
#include "utility/rng.h"
#include "utility/momentum.h"
#include "type_cast.h"

namespace Python {
bool Convert(Object obj, std::string& val)
{
    if (PyString_Check(obj.Get()))
        val = PyString_AsString(obj.Get());
    else if (PyUnicode_Check(obj.Get())) {
        Object sobj = PyUnicode_AsUTF8String(obj.Get());
        if (sobj.Get() == nullptr)
            return false;
        val = PyString_AsString(sobj.Get());
    }
    else
        return false;
    return true;
}
bool Convert(Object obj, bool& value)
{
    if (obj.Get() == Py_False)
        value = false;
    else if (obj.Get() == Py_True)
        value = true;
    else
        return false;
    return true;
}
bool IsInt(const Object& obj)
{
    return PyInt_Check(obj.Get()) || PyLong_Check(obj.Get());
}

bool Convert(Object obj, int& value)
{
    if (!IsInt(obj))
        return false;
    value = (int)PyLong_AsLong(obj.Get());
    return true;
}
bool Convert(Object obj, unsigned int& value)
{
    if (!IsInt(obj))
        return false;
    value = (unsigned int)PyLong_AsUnsignedLong(obj.Get());
    return true;
}
bool Convert(Object obj, long& value)
{
    if (!IsInt(obj))
        return false;
    value = PyLong_AsLong(obj.Get());
    return true;
}
bool Convert(Object obj, unsigned long& value)
{
    if (!IsInt(obj))
        return false;
    value = PyLong_AsUnsignedLong(obj.Get());
    return true;
}
bool Convert(Object obj, long long& value)
{
    if (!IsInt(obj))
        return false;
    value = PyLong_AsLongLong(obj.Get());
    return true;
}
bool Convert(Object obj, unsigned long long& value)
{
    if (!IsInt(obj))
        return false;
    value = PyLong_AsUnsignedLongLong(obj.Get());
    return true;
}

template <typename T>
bool ConvertReal(Object obj, T& val)
{
    if (!PyFloat_Check(obj.Get()) && !IsInt(obj))
        return false;
    val = (T)PyFloat_AsDouble(obj.Get());
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
    if (!PyComplex_Check(obj.Get()))
        return false;
    val = Complex(PyComplex_RealAsDouble(obj.Get()),
                  PyComplex_ImagAsDouble(obj.Get()));
    return true;
}
bool Convert(Object obj, RandomFactory& rng)
{
    Object sourceobj = PyObject_Str(obj.Get());
    std::string source;
    if (!Convert(sourceobj, source))
        return false;
    istringstream iss(source);
    iss >> rng;
    return !(iss.bad() || iss.fail());
}
bool Convert(Object obj, ITypeCast& value)
{
    return value.FromPy(obj);
}

bool Convert(Object obj, spin& val)
{
    return Convert(obj, val);
}
bool Convert(Object obj, Momentum& val)
{
    return Convert(obj, val.K);
}
// Allocation methods

Object CastToPy(const std::string& str)
{
    return PyString_FromString(str.c_str());
}
Object CastToPy(const char* str)
{
    return PyString_FromString(str);
}
Object CastToPy(bool value)
{
    return PyBool_FromLong(value);
}
Object CastToPy(int num)
{
    return PyInt_FromLong(num);
}
Object CastToPy(unsigned int num)
{
    return PyInt_FromLong(num);
}
Object CastToPy(long num)
{
    return PyLong_FromLong(num);
}
Object CastToPy(unsigned long num)
{
    return PyLong_FromUnsignedLong(num);
}
Object CastToPy(unsigned long long num)
{
    return PyLong_FromUnsignedLongLong(num);
}
Object CastToPy(long long num)
{
    return PyLong_FromLongLong(num);
}

Object CastToPy(float num)
{
    return PyFloat_FromDouble(num);
}
Object CastToPy(double num)
{
    return PyFloat_FromDouble(num);
}
Object CastToPy(const Complex& num)
{
    return PyComplex_FromDoubles(num.Re, num.Im);
}
Object CastToPy(const RandomFactory& rng)
{
    return CastToPy(ToString(rng));
}
Object CastToPy(const ITypeCast& val)
{
    return val.ToPy();
}

Object CastToPy(Object obj)
{
    return obj.Copy();
}

Object CastToPy(spin num)
{
    return CastToPy((int)num);
}
Object CastToPy(Momentum k)
{
    return CastToPy(k.K);
}
}
