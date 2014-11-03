//
//  parameter_map.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/1/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "parameter_map.h"
#include <iostream>
#include "../utility/scopeguard.h"

using namespace std;

#define SEP "#"
/**
*  Load all __existed__ keys from file
*
*  @param InputFile Input file name
*/
void ParameterMap::LoadFromFile(const std::string &InputFile)
{
    clear();
    ifstream ifs(InputFile, ios::in);
    ON_SCOPE_EXIT([&] {ifs.close(); });
    if (!ifs.is_open())
        ABORT("Fail to open input file " << InputFile);
    string temp;
    _map.clear();
    while (getline(ifs, temp)) {
        _map.insert(make_pair(temp));
    }
}
/**
*  Save all key/value pair to Outputfile
*
*  @param OutputFile OutputFile
*  @param Mode       "w" to write/"a" to append
*
*/
void ParameterMap::SaveToFile(const std::string &OutputFile, string Mode)
{
    auto mode = ios::out;
    if ((Mode) == "a")
        mode = ios::app;
    else if ((Mode) != "w")
        ABORT("I don't know what is the mode " << Mode << "?");

    ofstream ofs(OutputFile, mode);
    ON_SCOPE_EXIT([&] {ofs.close(); });
    if (!ofs.is_open())
        ABORT("Fail to open file " << OutputFile);
    for (auto &kv : _map)
        if (kv.second != "")
            ofs << kv.second << "    #" << kv.first << std::endl;
}

void ParameterMap::clear()
{
    _map.clear();
}

void ParameterMap::addKey(std::string key)
{
    _map.at(_ToUpper(key)) = "";
}

void ParameterMap::eraseKey(std::string key)
{
    if (_DoesKeyExist(key))
        _map.erase(_ToUpper(key));
}

std::string ParameterMap::get(std::string key)
{
    _MakeSureKeyExists(key);
    return _map[_ToUpper(key)];
}

string trim(string s)
{
    if (s.empty()) {
        return s;
    }

    s.erase(0, s.find_first_not_of(" "));
    s.erase(s.find_last_not_of(" ") + 1);
    return s;
}

std::pair<string, string> ParameterMap::make_pair(string source)
{
    auto pos = source.find(SEP);
    if (pos == string::npos)
        ABORT("Are you sure the separator is " << SEP << "?");
    string key = trim(source.substr(pos + 1, source.size()));
    return pair<string, string>(_ToUpper(key),
                                trim(source.substr(0, pos)));
}
string ParameterMap::_ToUpper(string source)
{
    for (auto &c : source)
        c = toupper(c);
    return source;
}

bool ParameterMap::_DoesKeyExist(string key)
{
    key = _ToUpper(key);
    if (_map.find(key) == _map.end())
        return false;
    else
        return true;
}
bool ParameterMap::_MakeSureKeyExists(string key)
{
    key = _ToUpper(key);
    if (_map.find(key) == _map.end()) {
        ABORT("Can not find the key " << key);
        return false;
    }
    return true;
}