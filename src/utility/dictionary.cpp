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

Object Dictionary::CastToPy() const
{
    return GetObject().Copy();
}
bool Dictionary::Convert(Object obj)
{
    if (!PyDict_Check(obj.Get()))
        return false;
    _Dict = obj;
    return true;
}
Dictionary::Dictionary()
{
    _Dict = PyDict_New();
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
    SaveDict.CallFunction("SaveDict", FileName, Mode, _Dict);
}
void Dictionary::Clear()
{
    PyDict_Clear(_Dict.Get());
}

void Dictionary::Print()
{
    _Dict.Print();
}

void Dictionary::_PrintDebug() const
{
    _Dict._PrintDebug();
}