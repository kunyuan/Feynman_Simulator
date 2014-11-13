//
//  control.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/3/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "message.h"

bool para::Message::Load()
{
    //do not abort if message file does not exist
    if (_Para.ParseFile(_MessageFile, false)) {
        GetParaArray(_Para, Interaction, MODEL_PARA_NUM);
        GetPara(_Para, ExternalField);
        GetPara(_Para, Beta);
        GetPara(_Para, Version);
        return true;
    }
    else
        return false;
}

void para::Message::Save()
{
    _Para.clear();
    SetParaArray(_Para, Interaction, MODEL_PARA_NUM);
    SetPara(_Para, ExternalField);
    SetPara(_Para, Beta);
    SetPara(_Para, Version);
    _Para.SaveToFile(_MessageFile, "w");
}

std::string para::Message::PrettyString()
{
    return _Para.PrettyString();
}
