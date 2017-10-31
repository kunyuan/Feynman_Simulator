//
//  control.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/3/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "message.h"
#include "utility/dictionary.h"

bool para::Message::Load(const string& FileName)
{
    //do not abort if message file does not exist
    Dictionary _Para;
    try {
        _Para.Load(FileName);
    }
    catch (IOInvalid e) {
        return false;
    }
    GET(_Para, Beta);
    GET(_Para, Version);
    GET(_Para, SqueezeFactor);
    return true;
}

void para::Message::Save(const string& FileName)
{
    Dictionary _Para;
    SET(_Para, Beta);
    SET(_Para, Version);
    SET(_Para, SqueezeFactor);
    _Para.Save(FileName, "w");
}

std::string para::Message::PrettyString()
{
    Dictionary _Para;
    SET(_Para, Beta);
    SET(_Para, Version);
    SET(_Para, SqueezeFactor);
    return _Para.PrettyString();
}
