//
//  main.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/2/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

/********************** define the running mode here **************************/
#define DEBUGLEVEL 0

/********************** include files *****************************************/
#include <iostream>
#include "initialization.h"
#include "definition_global.h"
using namespace std;
int counter;

int main(int argc, const char *argv[])
{
    LOGGER_CONF("", "MC", Logger::file_on | Logger::screen_on, INFO, INFO);

    Initilization();
    counter = 0;

    TestAll();
    return 0;
}
