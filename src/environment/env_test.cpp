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
    EnvMonteCarlo env(0);
    string inputfile = "../src/module/parameter/_in_MC_test";
    env.BuildNew(inputfile, true);
    sput_fail_unless(Equal(env.Para.OrderReWeight[3], 4.0), "check reading the job file");

    //do some change on env
    env.Para.RNG.urn();
    env.Para.OrderReWeight[3] = 4.1;
    //then save the state
    env.Save();

    EnvMonteCarlo new_env(0);
    new_env.BuildNew(inputfile, true);
    new_env.Load();
    sput_fail_unless(Equal(new_env.Para.OrderReWeight[3], 4.1), "Check reading state file");
    sput_fail_unless(Equal(env.Para.RNG.urn(), new_env.Para.RNG.urn()), "Check reading RNG from state file");
    new_env.DeleteSavedFiles();
}