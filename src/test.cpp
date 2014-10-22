//
//  test.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/4/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "test.h"
#include "abort.h"
using namespace std;

#define TEST(func)                  \
    {                               \
        if (func() != EXIT_SUCCESS) \
            ABORT("Shit happens!"); \
    }
int RunTest()
{

    //    TestTimer();  //Test the timer
    //        TestRNG();
    TEST(TestDiagram);
    //    TEST(TestLattice);
    //    TEST(TestObservable);

    //    TestArray();
    //    Testcnpy();
    return 0;
}