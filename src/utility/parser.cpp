//
//  parameter_map.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/1/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "parser.h"
#include <iostream>
#include "utility/scopeguard.h"

using namespace std;

#define SEP '#'

string trim(string s)
{
    if (s.empty()) {
        return s;
    }

    s.erase(0, s.find_first_not_of(" "));
    s.erase(s.find_last_not_of(" ") + 1);
    return s;
}

/**
*  Load all __existed__ keys from file
*
*  @param InputFile Input file name
*/
bool SimpleParser::ParseFile(const std::string &InputFile, bool AbortIfFail)
{
    clear();
    ifstream ifs(InputFile, ios::in);
    ON_SCOPE_EXIT([&] {ifs.close(); });
    if (!ifs.is_open()) {
        if (AbortIfFail) {
            ABORT("Fail to open input file " << InputFile);
        }
        else {
            LOG_WARNING(InputFile << " does not exist!");
            return false;
        }
    }
    string temp;
    _map.clear();
    while (getline(ifs, temp)) {
        string key = trim(temp);
        if (temp[0] == SEP && temp[1] == SEP)
            continue;
        _map.insert(make_pair(temp));
    }
    return true;
}
/**
*  Save all key/value pair to Outputfile
*
*  @param OutputFile OutputFile
*  @param Mode       "w" to write/"a" to append
*
*/
void SimpleParser::SaveToFile(const std::string &OutputFile, string Mode)
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

string SimpleParser::PrettyString()
{
    stringstream ss;
    for (auto &kv : _map)
        if (kv.second != "")
            ss << kv.first << "=" << kv.second << std::endl;
    return ss.str();
}

void SimpleParser::clear()
{
    _map.clear();
}

void SimpleParser::addKey(std::string key)
{
    _map.at(_ToUpper(key)) = "";
}

void SimpleParser::eraseKey(std::string key)
{
    if (_DoesKeyExist(key))
        _map.erase(_ToUpper(key));
}

std::string SimpleParser::get(std::string key)
{
    _MakeSureKeyExists(key);
    return _map[_ToUpper(key)];
}

std::pair<string, string> SimpleParser::make_pair(string source)
{
    auto pos = source.find(SEP);
    if (pos == string::npos)
        ABORT("Are you sure the separator is " << SEP << "?");
    string key = trim(source.substr(pos + 1, source.size()));
    return pair<string, string>(_ToUpper(key),
                                trim(source.substr(0, pos)));
}
string SimpleParser::_ToUpper(string source)
{
    for (auto &c : source)
        c = toupper(c);
    return source;
}

bool SimpleParser::_DoesKeyExist(string key)
{
    key = _ToUpper(key);
    if (_map.find(key) == _map.end())
        return false;
    else
        return true;
}
bool SimpleParser::_MakeSureKeyExists(string key)
{
    key = _ToUpper(key);
    if (_map.find(key) == _map.end()) {
        ABORT("Can not find the key " << key);
        return false;
    }
    return true;
}