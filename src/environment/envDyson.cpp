//
//  enDyson.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/18/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "environment.h"

EnvDyson::EnvDyson(int pid)
    : Environment(pid)
{
}

bool EnvDyson::BuildNew(const string &InputFile)
{
    return true;
}

string EnvDyson::_FinalQuanFile()
{
    return "final_quantity.dat";
}
string EnvDyson::_FinalStatisFile()
{
    return "final_quantity_statis.npz";
}