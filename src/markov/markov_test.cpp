//
//  markov_test.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/20/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "../environment/environment.h"
#include "../utility/sput.h"
using namespace std;

void Test_CreateWorm();

int TestMarkov()
{
    sput_start_testing();
    sput_enter_suite("Test Updates:");
    sput_run_test(Test_CreateWorm);
    sput_finish_testing();
    return sput_get_return_value();
}

void Test_CreateWorm()
{
    EnvMonteCarlo Env(0);
    Env.BuildNew("../src/environment/_in_MC_test");
    LOG_INFO("Build Environment succeed!");
    Env.Diag.Load("../src/diagram/diagram_template.config");
    LOG_INFO("Load Diagram from config file succeed!");
    Markov &markov = Env.Grasshopper;
    int total = 0;
    for (int i = 0; i < 10; i++) {
        markov.CreateWorm();
        if (markov.Diag->Worm.Exist) {
            total += 1;
            markov.Diag->Save("diagram_template.config", "a");
        }
        markov.Diag->Load("../src/diagram/diagram_template.config");
    }
    LOG_INFO("Updates(Create Worm) are done!");
    sput_fail_unless(total == 10, "The accept ratio of CreateWorm");
}
