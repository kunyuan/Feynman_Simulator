//
//  main.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/2/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//


/********************** include files *****************************************/
#include <iostream>
#include "definition.h"
#include "initialization.h"
using namespace std;

void MonteCarlo(Jobs &Job);
void Dyson(Jobs& Job);
int main(int argc, const char *argv[])
{
    InitGlobalUtility();
    RunTest();
    
    Jobs Job;
    Job.Read();
    
    InitEveryOneNeedsIt();
    switch(Job.Type)
    {
        case MC:
            MonteCarlo(Job);
        case DYSON:
            Dyson(Job);
    }
    return 0;
}

void MonteCarlo(Jobs &Job)
{
    EnvMoneCarlo Env(Job);
}

void Dyson(Jobs &Job)
{
    EnvDyson Env(Job);
}
