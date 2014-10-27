//
//  markov_test.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/20/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include <sstream>
#include "markov.h"
#include "markov_monitor.h"
#include "sput.h"
#include "utility.h"
#include "diagram.h"
#include "logger.h"
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
    Jobs Job;
    EnvMonteCarlo Env(Job);
    Markov markov(&Env);
    markov.Diag->LoadConfig("../src/diagram/diagram_template.config");
    int total = 0;
    for (int i = 0; i < 10; i++) {
        markov.CreateWorm();
        if (markov.Diag->Worm.Exist) {
            total += 1;
            markov.Diag->SaveConfig("diagram_template.config", "w");
        }
        markov.Diag->LoadConfig("../src/diagram/diagram_template.config");
    }
    cout << total << endl;
    sput_fail_unless(total == 10, "The accept ratio of CreateWorm");
}
