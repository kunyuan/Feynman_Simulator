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

void Dictionary::LoadFromString(const std::string& script)
{
    AnyObject obj;
    obj.EvalScript(script);
    if (!FromPy(obj))
        ERRORCODEABORT(ERR_VALUE_INVALID, "Script is invalided!");
}

void Dictionary::Load(const std::string& FileName)
{
    if (!DoesFileExist(FileName))
        ERRORCODEABORT(ERR_FILE_NOT_FOUND, FileName + " does not exist!");
    ModuleObject LoadDict;
    LoadDict.LoadModule("IO.py");
    Object result = LoadDict.CallFunction("LoadDict", FileName);
    if (!FromPy(result))
        ERRORCODEABORT(ERR_VALUE_INVALID, "File is invalided!");
}

void Dictionary::Save(const string& FileName, const std::string& Mode,
                      const std::string& key)
{
    ModuleObject SaveDict;
    SaveDict.LoadModule("IO.py");
    if (Mode == "a") {
        Dictionary root;
        try {
            root.Load(FileName);
        }
        catch (ERRORCODE e) {
            if (e != ERR_FILE_NOT_FOUND)
                throw e;
        }
        if (root.IsEmpty())
            root = *this;
        else
            root[key] = *this;
        root.Save(FileName, "w");
    }
    else
        SaveDict.CallFunction("SaveDict", FileName, "w", _Map);
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