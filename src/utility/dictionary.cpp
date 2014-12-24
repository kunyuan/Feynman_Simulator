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

bool Dictionary::FromPy(Object obj)
{
    return Python::Convert(obj, _Map);
}

void Dictionary::LoadFromString(const std::string& script)
{
    AnyObject obj;
    obj.EvalScript(script);
    if (!FromPy(obj))
        ERRORCODEABORT(ERR_VALUE_INVALID, "Script is invalided!");
}

void Dictionary::Load(const std::string& FileName)
{
    ModuleObject LoadDict;
    LoadDict.LoadModule("IO.py");
    Object result = LoadDict.CallFunction("LoadDict", FileName);
    if (!FromPy(result))
        ERRORCODEABORT(ERR_VALUE_INVALID, "File is invalided!");
}

void Dictionary::Save(const string& FileName, const std::string& Mode)
{
    ModuleObject SaveDict;
    SaveDict.LoadModule("IO.py");
    SaveDict.CallFunction("SaveDict", FileName, "w", _Map);
}
void Dictionary::Clear()
{
    _Map.clear();
}

void Dictionary::Print()
{
    AnyObject(_Map).Print();
}

std::string Dictionary::PrettyString()
{
    return AnyObject(_Map).PrettyString();
}