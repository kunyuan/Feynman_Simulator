//
//  enDyson.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/18/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "environment.h"

EnvDyson::EnvDyson()
{
}

bool EnvDyson::BuildFromFile(string InputFile)
{
    Environment::BuildFromFile(InputFile);
    LOGGER_CONF(_LogFile(), "DYSON", Logger::file_on | Logger::screen_on, INFO, INFO);
    return true;
}

string EnvDyson::_TotalStatisFile()
{
    return "global_statis.npz";
}
string EnvDyson::_FinalQuanFile()
{
    return "final_quantity.dat";
}
string EnvDyson::_FinalStatisFile()
{
    return "final_quantity_statis.npz";
}