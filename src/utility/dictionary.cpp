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
#include "dictionary.h"

using namespace std;
using namespace Python;

Object Dictionary::ToPy() const
{
    return Copy();
}
bool Dictionary::FromPy(Object obj)
{
    if (!PyDict_Check(obj.Get()))
        return false;
    *this = obj;
    return true;
}

void Dictionary::LoadFromString(const std::string& script)
{
    AnyObject obj;
    obj.EvalScript(script);
    *this = obj;
}

void Dictionary::Load(const std::string& FileName)
{
    ModuleObject LoadDict;
    LoadDict.LoadModule("IO.py");
    Object result = LoadDict.CallFunction("LoadDict", FileName);
    *this = result;
}

void Dictionary::Save(const string& FileName, const string& Mode)
{
    ModuleObject SaveDict;
    SaveDict.LoadModule("IO.py");
    SaveDict.CallFunction("SaveDict", FileName, Mode, *this);
}
void Dictionary::Clear()
{
    PyDict_Clear(_PyPtr);
}