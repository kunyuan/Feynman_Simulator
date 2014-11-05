//
//  control.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/3/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "status.h"

bool para::status::Load()
{
    //do not abort if Status file does not exist
    if (_Para.ParseFile(_StatusFile, false)) {
        GetPara(_Para, Jcp);
        GetPara(_Para, Beta);
        GetPara(_Para, Version);
        GetPara(_Para, OrderAccepted);
        return true;
    }
    else
        return false;
}

void para::status::Save()
{
    _Para.clear();
    SetPara(_Para, Jcp);
    SetPara(_Para, Beta);
    SetPara(_Para, Version);
    SetPara(_Para, OrderAccepted);
    _Para.SaveToFile(_StatusFile, "w");
}

std::string para::status::PrettyString()
{
    return _Para.PrettyString();
}
