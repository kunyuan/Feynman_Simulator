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
bool Convert(Object obj, Dictionary& value)
{
    if (!PyDict_Check(obj.Borrow()))
        return false;
    value = Dictionary(obj);
    return true;
}
Object CastToPy(const Dictionary& value)
{
    //    PyObject* dptr = value.GetObject().Borrow();
    //    Py_INCREF(dptr);
    //    return Object::Steal(dptr);
    return value.GetObject();
}
}

Dictionary::Dictionary()
{
    _Dict = Object::Steal(PyDict_New());
}

void Dictionary::LoadFromString(const std::string& script)
{
    _Dict.EvalScript(script);
}

void Dictionary::Load(const std::string& FileName)
{
    ModuleObject LoadDict;
    LoadDict.LoadModule("IO.py");
    Object result = LoadDict.CallFunction("LoadDict", FileName);
    _Dict = result;
}

void Dictionary::Save(const string& FileName, const string& Mode)
{
    ModuleObject SaveDict;
    SaveDict.LoadModule("IO.py");
    //    SaveDict.CallFunction("Simple");
    SaveDict.CallFunction("SaveDict", FileName, Mode, _Dict);
}
void Dictionary::Clear()
{
    PyDict_Clear(_Dict.Borrow());
}

void Dictionary::Print()
{
    _Dict.Print();
}

void Dictionary::_PrintDebug() const
{
    LOG_INFO("PyObject ref=" << _Dict.Borrow()->ob_refcnt);
}