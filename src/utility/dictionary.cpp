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
    ASSERT_ALLWAYS(DoesFileExist(FileName), "File not exist!");
    ifstream ifs(FileName, std::ios::in);
    ON_SCOPE_EXIT([&] {ifs.close(); });
    string source;
    ifs >> source;
    LoadByEval(source);
}

void Dictionary::Save(const string& FileName, const string& Mode)
{
    auto mode = ios::out;
    if ((Mode) == "a")
        mode = ios::app;
    else if ((Mode) != "w")
        ABORT("I don't know what is the mode " << Mode << "?");

    ofstream ofs(FileName, mode);
    ON_SCOPE_EXIT([&] {ofs.close(); });
    if (!ofs.is_open())
        ABORT("Fail to open file " << FileName);
    ofs << _Dict.PrettyString() << endl;
}

void Dictionary::Print()
{
    _Dict.Print();
}

void Dictionary::_PrintDebug() const
{
    LOG_INFO("Object ref=" << _Dict.use_count() << ", PyObject ref=" << _Dict.get()->ob_refcnt);
}