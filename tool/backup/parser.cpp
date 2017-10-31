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

#define SEP '='
#define COMMENT '#'

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
bool SimpleParser::ParseFile(const std::string& InputFile)
{
    clear();
    ifstream ifs(InputFile, ios::in);
    ON_SCOPE_EXIT([&] {ifs.close(); });
    if (!ifs.is_open()) {
        LOG_WARNING(InputFile << " does not exist!");
        throw(ERR_FILE_NOT_FOUND);
        return false;
    }
    string temp;
    _map.clear();
    while (getline(ifs, temp)) {
        string key = trim(temp);
        if (temp.empty() || temp[0] == COMMENT)
            //empty line or comment line
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
void SimpleParser::SaveToFile(const std::string& OutputFile, string Mode)
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
    for (auto& kv : _map)
        if (kv.second != "")
            ofs << kv.first << " = " << kv.second << std::endl;
}

string SimpleParser::PrettyString()
{
    stringstream ss;
    for (auto& kv : _map)
        if (kv.second != "")
            ss << kv.first << "=" << kv.second << std::endl;
    return ss.str();
}

void SimpleParser::clear()
{
    _map.clear();
}

std::pair<string, string> SimpleParser::make_pair(string source)
{
    auto pos = source.find(SEP);
    if (pos == string::npos)
        ABORT("Are you sure the separator is " << SEP << "?");
    string key = trim(source.substr(0, pos));
    return pair<string, string>(key, trim(source.substr(pos + 1, source.size())));
}

bool SimpleParser::_DoesKeyExist(string key)
{
    if (_map.find(key) == _map.end())
        return false;
    else
        return true;
}
bool SimpleParser::_MakeSureKeyExists(string key)
{
    if (_map.find(key) == _map.end()) {
        ABORT("Can not find the key " << key);
        return false;
    }
    return true;
}