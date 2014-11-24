//
//  weight0_test.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/23/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "component.h"
#include "utility/sput.h"
#include "utility/rng.h"

using namespace std;
using namespace weight0;

void TestWeightMatrix();
int weight0::TestWeight()
{
    sput_start_testing();
    sput_enter_suite("Test Weight...");

    //test diagram object weight, like sigma, G
    sput_run_test(TestWeightMatrix);
    sput_finish_testing();
    return sput_get_return_value();
}

void TestWeightMatrix()
{
    Lattice lat({8, 8}, CHECKBOARD);
    real Beta = 1.0;
    G G1(J1J2, lat, Beta);
    G1.BuildNew();
    G G2 = G1;
    G2.Save("test_weight.npz");
}