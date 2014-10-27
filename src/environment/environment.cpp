//
//  environment.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/17/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "environment.h"

Environment::Environment(const JobsBase &Job)
    : Lat(), Beta(Job.Beta), Order(Job.Order)
{
    StateFile = "test.dat";
}

void Environment::SaveState()
{
}
