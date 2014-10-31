//
//  env_test.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/20/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "environment.h"
#include "../utility/sput.h"
#include "../utility/utility.h"

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
    EnvMonteCarlo env;
    string inputfile = "../src/environment/_in_MC_test";
    env.BuildFromFile(inputfile);
    sput_fail_unless(Equal(env.OrderWeight[3], 4.0), "check reading the job file");

    //do some change on env
    env.RNG.urn();
    env.OrderWeight[3] = 4.1;
    //then save the state
    env.SaveState();

    EnvMonteCarlo new_env;
    new_env.PID = env.PID;
    new_env.LoadState();
    sput_fail_unless(Equal(new_env.OrderWeight[3], 4.1), "Check reading state file");
    sput_fail_unless(Equal(env.RNG.urn(), new_env.RNG.urn()), "Check reading RNG from state file");
}