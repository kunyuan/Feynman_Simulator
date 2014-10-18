//
//  initialization.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/3/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "initialization.h"
using namespace std;

void InitGlobalUtility()
{
    //initialize LOGGER
    LOGGER_CONF("", "MC", Logger::file_on | Logger::screen_on, INFO, INFO);
    //initialize Counter
    Counter=0;
    //Randomize
    RNG.Reset();
}

void InitEveryOneNeedsIt()
{
    
}

void InitMonteCarlo()
{
    return;
}

void InitDyson()
{
    return;
}
