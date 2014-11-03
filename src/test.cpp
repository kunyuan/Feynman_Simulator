//
//  test.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/4/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "test.h"
#include "utility/crc32.h"
#include "markov/markov.h"
#include "diagram/diagram.h"
#include "lattice/lattice.h"
#include "observable/weight.h"
#include "utility/fft.h"
#include "environment/environment.h"
using namespace std;

#define TEST(func)                  \
    {                               \
        if (EXIT_SUCCESS != func()) \
            exit(0);                \
    }
int RunTest()
{

    //    TestTimer();  //Test the timer
    //    TestRNG();
    //    TestArray();
    //    TEST(TestEnvironment);
    //    TEST(TestDiagram);
    //    TEST(TestLattice);
    //    TEST(TestObservable);
    TEST(TestMarkov);

    //    TEST(TestCRC32);
    //    TEST(TestFFT);

    //    Testcnpy();
    return 0;
}
