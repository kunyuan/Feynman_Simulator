//
//  control.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/3/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "message.h"
#include "utility/dictionary.h"

bool para::Message::Load()
{
    //do not abort if message file does not exist
    Dictionary _Para;
    try {
        _Para.Load(_MessageFile);
    }
    catch (ERRORCODE e) {
        if (e != ERR_FILE_NOT_FOUND)
            throw e;
        return false;
    }
    GET(_Para, Beta);
    GET(_Para, Version);
    return true;
}

void para::Message::Save()
{
    Dictionary _Para;
    SET(_Para, Beta);
    SET(_Para, Version);
    _Para.Save(_MessageFile, "w");
}

std::string para::Message::PrettyString()
{
    Dictionary _Para;
    SET(_Para, Beta);
    SET(_Para, Version);
    return _Para.PrettyString();
}
