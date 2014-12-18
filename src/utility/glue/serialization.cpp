//
//  serialization.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 12/17/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "serialization.h"
#include "utility/complex.h"
#include <vector>
#include <string>
#include "glue.h"
#include <initializer_list>
#include <iostream>

using namespace std;
Serialization::Serialization()
{
    _Dict = NewDict();
    LOG_INFO("reference cout=" << _Dict->ob_refcnt);
}

void Serialization::Set(std::string key, int value)
{
    //    auto object = PyInt_FromLong(value);
    //    LOG_INFO("reference cout=" << object->ob_refcnt);
    //    PyDict_SetItemString(_Dict, key.c_str(), object);
    //    LOG_INFO("reference cout=" << object->ob_refcnt);
    //    Py_DECREF(object);
    DictSet(_Dict, key, CastInt(value));
}

int Serialization::Get(std::string key)
{
    //    auto object = PyDict_GetItemString(_Dict, key.c_str());
    //    return (int)PyInt_AsLong(object);
    int result;
    AsInt(DictGet(_Dict, key), result);
    return result;
}

void Serialization::Print()
{
    DictPrint(_Dict);
    //    PyObject_Print(_Dict, stdout, 0);
    //    PyObject* main = PyImport_AddModule("__main__");
    //    PyObject* globalDictionary = PyModule_GetDict(main);
    //    PyObject* localDictionary = PyDict_New();
    //    const char* pythonScript = "{'a':1,'b':2}";
    //    //    PyDict_SetItemString(localDictionary, "multiplicand", PyInt_FromLong(2));
    //    //    PyDict_SetItemString(localDictionary, "multiplier", PyInt_FromLong(5));
    //    PyObject* result = PyRun_String(pythonScript, Py_eval_input, globalDictionary, localDictionary);
    //    PyObject_Print(result, stdout, 0);
}

Serialization::~Serialization()
{
    Py_DECREF(_Dict);
    LOG_INFO("reference cout=" << _Dict->ob_refcnt);
}