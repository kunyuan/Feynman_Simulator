//
//  env_test.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/20/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "environment.h"
#include "utility/sput.h"

void Test_EnvMC();

int TestEnvironment()
{
    sput_start_testing();
    sput_enter_suite("Test Environment");
    sput_run_test(Test_EnvMC);
    sput_finish_testing();
    return sput_get_return_value();
}

void Test_EnvMC()
{
}