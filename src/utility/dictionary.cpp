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
    return Python::CastToPy(_Map);
}

bool Dictionary::FromPy(const Object& obj)
{
    return Python::Convert(obj, _Map);
}

Python::AnyObject& Dictionary::operator[](const std::string& key)
{
    return _Map[key];
}

void Dictionary::Update(const Dictionary& dict)
{
    for (auto e : dict)
        _Map[e.first] = e.second;
}
void Dictionary::LoadFromString(const std::string& script)
{
    AnyObject obj;
    obj.EvalScript(script);
    if (!FromPy(obj))
        ABORT("Script is invalided!");
}

void Dictionary::Load(const std::string& FileName)
{
    ModuleObject LoadDict;

    LoadDict.LoadModule("IO.py");
    Object result = LoadDict.CallFunction("LoadDict", FileName);
    PropagatePyError();
    if (!FromPy(result))
        ABORT("Fail to read file!");
}

void Dictionary::Save(const string& FileName, const std::string& Mode)
{
    ModuleObject SaveDict;
    SaveDict.LoadModule("IO.py");
    SaveDict.CallFunction("SaveDict", FileName, Mode, _Map);
    PropagatePyError();
}

void Dictionary::BigLoad(const std::string& FileName)
{
    ModuleObject LoadBigDict;
    LoadBigDict.LoadModule("IO.py");
    Object result = LoadBigDict.CallFunction("LoadBigDict", FileName);
    PropagatePyError();
    if (!FromPy(result))
        ABORT("Fail to read file " << FileName);
}
void Dictionary::BigSave(const std::string& FileName)
{
    ModuleObject SaveBigDict;
    SaveBigDict.LoadModule("IO.py");
    SaveBigDict.CallFunction("SaveBigDict", FileName, _Map);
    PropagatePyError();
}

void Dictionary::Clear()
{
    _Map.clear();
}

bool Dictionary::HasKey(const std::string& key) const
{
    return _Map.find(key) != _Map.end();
}
bool Dictionary::IsEmpty() const
{
    return _Map.empty();
}
void Dictionary::Print() const
{
    AnyObject(_Map).Print();
}

std::string Dictionary::PrettyString() const
{
    return AnyObject(_Map).PrettyString();
}