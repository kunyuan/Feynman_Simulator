//
//  dyson_test.cpp
//  Feynman_Simulator
//
//  Created by yuan on 11/19/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include <stdio.h>
#include "dyson.h"
#include "module/weight/weight.h"
#include "module/weight/weight_inherit.h"
#include "module/parameter/parameter.h"
#include "environment/environment.h"
#include "utility/sput.h"
#include "utility/utility.h"
#include "utility/rng.h"

using namespace std;
using namespace weight;
using namespace dyson;

void TestMultiply();
void TestInverse();
void TestG();
void TestW();

int dyson::TestDyson()
{
    sput_start_testing();
    sput_enter_suite("Test Dyson...");
    
    sput_run_test(TestMultiply);
    sput_run_test(TestInverse);
    sput_run_test(TestG);
    sput_run_test(TestW);

    sput_finish_testing();
    return sput_get_return_value();
}

void TestMultiply()
{
    para::ParaDyson Para;
    Para.SetTest();
    Weight Weight(true);
    Weight.SetTest(Para);
    
    unsigned int GShape[5];
    AssignFromTo(&GShape[SP], Weight.G->Shape(), 4);
    
    unsigned int WShape[5];
    AssignFromTo(&WShape[SP], Weight.W->Shape(), 4);
    
    sput_fail_unless(GShape[SP]==4, "Check: dimension of spin of G");
    sput_fail_unless(WShape[SP]==16, "Check: dimension of spin of W");
    sput_fail_unless(GShape[SUB]==NSublattice2, "Check: dimension of sublattice of G");
    sput_fail_unless(WShape[SUB]==NSublattice2, "Check: dimension of sublattice of W");
    
    int spin_up = Weight.G->SpinIndex(UP, UP);
    int sub = Para.Lat.Sublat2Index(0, 0);
    int corr = Para.Lat.Vec2Index({0,0});
    int randomt = 8;
    Complex oldvalue = Weight.G->SmoothWeight[spin_up][sub][corr][randomt];
    MatrixInverse(Weight.G->SmoothWeight[spin_up], GShape[VOL]*GShape[TAU]);
    sput_fail_unless(Equal(Weight.G->SmoothWeight[spin_up][sub][corr][randomt], 1.0/oldvalue),
                     "Check: matrix inverse of G");
}

void TestInverse()
{
    
}

void TestG()
{
    
}

void TestW()
{
    
}
