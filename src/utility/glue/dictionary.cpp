//
//  serialization.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 12/17/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "dictionary.h"
#include "utility/complex.h"
#include <stdio.h>
#include "glue.h"
#include <initializer_list>
#include <iostream>
#include <fstream>
#include <memory>
#include "utility/abort.h"

using namespace std;

Dictionary::Dictionary()
{
    _Dict = PyDict_New();
    LOG_INFO("reference cout=" << _Dict->ob_refcnt);
}

Dictionary::~Dictionary()
{
    Py_XDECREF(_Dict);
    LOG_INFO("reference cout=" << _Dict->ob_refcnt);
}

bool Dictionary::LoadByEval(const std::string& source)
{
    PyObject* strobject = PyString_FromString(source.c_str());
    ON_SCOPE_EXIT([&] {Py_XDECREF(strobject); });
    PyObject* main = PyImport_AddModule("__main__");
    ON_SCOPE_EXIT([&] {Py_XDECREF(main); });
    PyObject* global = PyModule_GetDict(main); //borrowed reference
    PyObject* local = PyDict_New();
    ON_SCOPE_EXIT([&] {Py_XDECREF(local); });
    PyObject* dict = PyRun_String(source.c_str(), Py_eval_input,
                                  global, local);
    _CheckError();
    ASSERT_ALLWAYS(PyDict_Check(dict),
                   "a python dictionary is expected from the input source string!");
    Py_XDECREF(_Dict);
    _Dict = dict;
    return true;
}

bool Dictionary::Load(const std::string& FileName)
{
    ASSERT_ALLWAYS(DoesFileExist(FileName), "File not exist!");
    ifstream ifs(FileName, std::ios::in);
    ON_SCOPE_EXIT([&] {ifs.close(); });
    string source;
    ifs >> source;
    return LoadByEval(source);
}

void Dictionary::Save(const string& FileName, const string& Mode)
{
    FILE* fp = fopen(FileName.c_str(), Mode.c_str());
    ON_SCOPE_EXIT([&] {fclose(fp); });
    ASSERT_ALLWAYS(fp != nullptr, "File " << FileName << " fail to write!");
    PyObject_Print(_Dict, fp, 0);
}

void Dictionary::Print()
{
    PyObject_Print(_Dict, stdout, 0);
}

template <>
PyObject* Dictionary::_Cast(const Dictionary& value)
{
    return _Dict;
}

template <>
PyObject* Dictionary::_Cast(const int& value)
{
    PyObject* object = PyInt_FromLong(value);
    CHECKNULL(object);
    return object;
}

template <>
void Dictionary::_GetAs(PyObject* pyvalue, int& value)
{
    CHECKNULL(pyvalue);
    auto temp = PyInt_AsLong(pyvalue);
    _CheckError();
    value = (int)temp;
}

void Dictionary::_CheckError()
{
    if (PyErr_Occurred()) {
        PyErr_Print();
        ABORT("Python fatal error!");
    }
}