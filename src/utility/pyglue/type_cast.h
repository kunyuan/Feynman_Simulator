//
//  type_cast.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 12/20/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__type_cast__
#define __Feynman_Simulator__type_cast__

#include <string>
#include <utility>
#include <sstream>
#include <memory>
#include <map>
#include <vector>
#include <list>
#include <tuple>
#include "utility/vector.h"
#include "type_cast_interface.h"
#include "utility/convention.h"
#include <Python.h>

class Complex;
class RandomFactory;
class Momentum;
namespace Python {
// ------------ Conversion functions ------------

// Convert a PyObject to a std::string.
bool Convert(Object obj, std::string& val);
bool Convert(Object obj, char* val);
// Convert a PyObject to a bool value.
bool Convert(Object obj, bool& value);
bool Convert(Object obj, int& value);
bool Convert(Object obj, unsigned int& value);
bool Convert(Object obj, long& value);
bool Convert(Object obj, unsigned long& value);
bool Convert(Object obj, long long& value);
bool Convert(Object obj, unsigned long long& value);
bool Convert(Object obj, float& value);
bool Convert(Object obj, double& value);
bool Convert(Object obj, Complex& val);
bool Convert(Object obj, RandomFactory& val);
bool Convert(Object obj, spin& val);
bool Convert(Object obj, Momentum& val);
bool Convert(Object obj, ITypeCast& val);

template <typename T>
bool Convert(Object obj, Vec<T>& val)
{
    auto Dim = val.size();
    if (!PyList_Check(obj.Get()) || PyList_Size(obj.Get()) < Dim)
        return false;

    for (auto i = 0; i < Dim; i++) {
        T v;
        Object item = Object(PyList_GetItem(obj.Get(), i), NoRef);
        if (!Convert(item, v))
            return false;
        val[i] = v;
    }
    return true;
}
template <size_t n, class... Args>
typename std::enable_if<n == 0, bool>::type
AddToTuple(Object obj, std::tuple<Args...>& tup)
{
    Object item = Object(PyTuple_GetItem(obj.Get(), n), NoRef);
    return Convert(item, std::get<n>(tup));
}

template <size_t n, class... Args>
typename std::enable_if<n != 0, bool>::type
AddToTuple(Object obj, std::tuple<Args...>& tup)
{
    AddToTuple<n - 1, Args...>(obj.Get(), tup);
    Object item = Object(PyTuple_GetItem(obj.Get(), n), NoRef);
    return Convert(item, std::get<n>(tup));
}

template <class... Args>
bool Convert(Object obj, std::tuple<Args...>& tup)
{
    if (!PyTuple_Check(obj.Get()) || PyTuple_Size(obj.Get()) != sizeof...(Args))
        return false;
    return AddToTuple<sizeof...(Args)-1, Args...>(obj, tup);
}
// Convert a PyObject to a std::map
template <class K, class V>
bool Convert(Object obj, std::map<K, V>& mp)
{
    if (!PyDict_Check(obj.Get()))
        return false;
    PyObject* py_key, *py_val;
    Py_ssize_t pos(0);
    while (PyDict_Next(obj.Get(), &pos, &py_key, &py_val)) {
        //PyDict_Next return borrowed key and val
        K key;
        if (!Convert(Object(py_key, NoRef), key))
            return false;
        V val;
        if (!Convert(Object(py_val, NoRef), val))
            return false;
        mp.insert(std::make_pair(key, val));
    }
    return true;
}
// Convert a PyObject to a generic container.
template <class T, class C>
bool ConvertList(Object obj, C& container)
{
    if (!PyList_Check(obj.Get()))
        return false;
    for (Py_ssize_t i(0); i < PyList_Size(obj.Get()); ++i) {
        T val;
        Object item = Object(PyList_GetItem(obj.Get(), i), NoRef);
        if (!Convert(item, val))
            return false;
        container.push_back(std::move(val));
    }
    return true;
}
// Convert a PyObject to a std::list.
template <class T>
bool Convert(Object obj, std::list<T>& lst)
{
    return ConvertList<T, std::list<T> >(obj, lst);
}
// Convert a PyObject to a std::vector.
template <class T>
bool Convert(Object obj, std::vector<T>& vec)
{
    return ConvertList<T, std::vector<T> >(obj, vec);
}

// -------------- PyObject allocators ----------------

// Creates a PyObject from a std::string
Object CastToPy(const std::string& str);
Object CastToPy(const char* str);
Object CastToPy(int num);
Object CastToPy(unsigned int num);
Object CastToPy(long num);
Object CastToPy(unsigned long num);
Object CastToPy(long long num);
Object CastToPy(unsigned long long num);
// Creates a PyObject from a bool
Object CastToPy(bool value);
// Creates a PyObject from a real
Object CastToPy(float num);
Object CastToPy(double num);
Object CastToPy(const Complex& num);
Object CastToPy(const RandomFactory& val);
Object CastToPy(spin num);
Object CastToPy(Momentum num);
Object CastToPy(const ITypeCast& val);
Object CastToPy(Object);
// Creates a PyObject from a std::vector

// Generic python list allocation
template <class T>
static Object CastToPyList(const T& container)
{
    Object lst = PyList_New(container.size());

    Py_ssize_t i(0);
    for (auto it(container.begin()); it != container.end(); ++it) {
        auto obj = CastToPy(*it);
        PyList_SetItem(lst.Get(), i++, obj.Get(NewRef));
    }
    return lst;
}
template <typename T>
Object CastToPy(const Vec<T>& container)
{
    return CastToPyList(container);
}

template <class T>
Object CastToPy(const std::vector<T>& container)
{
    return CastToPyList(container);
}
// Creates a PyObject from a std::list
template <class T>
Object CastToPy(const std::list<T>& container)
{
    return CastToPyList(container);
}
// Creates a PyObject from a std::map
template <class T, class K>
Object CastToPy(
    const std::map<T, K>& container)
{
    Object dict = PyDict_New();

    for (auto it(container.begin()); it != container.end(); ++it) {
        Object first = CastToPy(it->first);
        Object second = CastToPy(it->second);
        PyDict_SetItem(dict.Get(), first.Get(), second.Get());
    }
    return dict;
}
}
#endif /* defined(__Feynman_Simulator__type_cast__) */
