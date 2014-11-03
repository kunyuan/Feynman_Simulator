//
//  control.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/3/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "state.h"

bool State::Load()
{
    if (_Para.LoadFromFile(_StateFile)) {
        GetPara(_Para, Jcp);
        GetPara(_Para, Beta);
        GetPara(_Para, Version);
        GetPara(_Para, OrderAccepted);
        return true;
    }
    else
        return false;
}

void State::Save()
{
    _Para.clear();
    SetPara(_Para, Jcp);
    SetPara(_Para, Beta);
    SetPara(_Para, Version);
    SetPara(_Para, OrderAccepted);
    _Para.SaveToFile(_StateFile, "w");
}
