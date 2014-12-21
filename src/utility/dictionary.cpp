//
//  serialization.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 12/17/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "utility/complex.h"
#include <initializer_list>
#include <fstream>
#include <memory>
#include "utility/abort.h"
#include "utility/scopeguard.h"
#include <iostream>
#include "dictionary.h"

using namespace std;
using namespace Python;

namespace Python {
bool Convert(PyObject* obj, Dictionary& value)
{
    if (!PyDict_Check(obj))
        return false;
    value = Dictionary(obj);
    return true;
}

PyObject* CastToPyObject(const Dictionary& value)
{
    PyObject* NewDict = value.GetPyObject();
    Py_INCREF(NewDict);
    return NewDict;
}
}

Dictionary::Dictionary()
{
    _Dict = PyDict_New();
}

void Dictionary::LoadByEval(const std::string& script)
{
    _Dict.EvalScript(script);
}

void Dictionary::Load(const std::string& FileName)
{
    Python::Object LoadDict;
    LoadDict.LoadScript("IO.py");
    auto result = LoadDict.CallFunction("LoadDict", FileName);
    _Dict = result;
}

void Dictionary::Save(const string& FileName, const string& Mode)
{
    Python::Object SaveDict;
    SaveDict.LoadScript("IO.py");
    auto result = SaveDict.CallFunction("SaveDict", FileName, Mode, _Dict);
}
void Dictionary::Clear()
{
    PyDict_Clear(_Dict.get());
}

void Dictionary::Print()
{
    _Dict.Print();
}

void Dictionary::_PrintDebug() const
{
    LOG_INFO("Object ref=" << _Dict.use_count() << ", PyObject ref=" << _Dict.get()->ob_refcnt);
}