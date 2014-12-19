//
//  serialization.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 12/17/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__serialization__
#define __Feynman_Simulator__serialization__

#include <Python/Python.h>
//#include "numpy/arrayobject.h"

#include <string>
#include <vector>
#include "utility/convention.h"
#include "utility/scopeguard.h"
#include "utility/utility.h"

class Dictionary {
public:
    Dictionary();
    ~Dictionary();
    template <typename T>
    void Set(const std::string& key, const T& value)
    {
        auto object = _Cast(value);
        ON_SCOPE_EXIT([&] {Py_XDECREF(object); });
        PyDict_SetItemString(_Dict, key.c_str(), object);
    }
    template <typename T>
    void Get(const std::string& key, T& value)
    {
        PyObject* object = PyDict_GetItemString(_Dict, key.c_str());
        _GetAs(object, value);
    }
    template <typename T>
    T Get(std::string key)
    {
        T value;
        Get(key, value);
        return value;
    }
    void Print();
    bool LoadByEval(const std::string&);
    bool Load(const std::string& FileName);
    void Save(const std::string& FileName, const std::string& Mode = "a");

private:
    PyObject* _Dict;
    void _CheckError();
    template <typename T>
    PyObject* _Cast(const T& value);
    template <typename T>
    PyObject* _Cast(const std::vector<T>& value)
    {
        PyObject* list = PyList_New(0);
        for (auto& e : value) {
            PyObject* temp = _Cast(e);
            ON_SCOPE_EXIT([&] {Py_XDECREF(temp); });
            PyList_Append(list, temp); //does not steal reference of PyObject* temp
        }
        return list;
    }
    template <typename T>
    void _GetAs(PyObject* pyvalue, T& value);
};

int TestDictionary();

#endif /* defined(__Feynman_Simulator__serialization__) */
