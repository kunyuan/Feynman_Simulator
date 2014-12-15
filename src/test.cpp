//
//  test.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/4/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "test.h"
#include "utility/crc32.h"
#include "environment/environment.h"
#include "module/markov/markov.h"
#include "module/diagram/diagram.h"
#include "module/weight/weight.h"
#include "lattice/lattice.h"
#include "estimator/estimator.h"
#include "module/weight/component.h"

using namespace std;

#define TEST(func)                  \
    {                               \
        if (EXIT_SUCCESS != func()) \
            exit(0);                \
    }
int RunTest()
{
    LOGGER_CONF("test.log", "test", Logger::file_on | Logger::screen_on, INFO, INFO);
    //    TestTimer();  //Test the timer
    //    TestRNG();
    //    TestArray();
    TEST(TestEnvironment);
    TEST(diag::TestDiagram);
    TEST(weight::TestWeight);
    TEST(mc::TestMarkov);
    TEST(mc::TestDiagCounter);

    //    TEST(TestLattice);
    //    TEST(TestEstimator);

    //    TEST(TestCRC32);
    //    TEST(TestFFT);

    //    Testcnpy();
    return 0;
}
