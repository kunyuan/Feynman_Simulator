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
    return true;
}