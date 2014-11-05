//
//  markov_test.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/20/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "markov.h"
#include "utility/sput.h"
#include "module/diagram/diagram.h"
#include "module/observable/weight.h"
#include "module/parameter/parameter.h"
using namespace std;
using namespace mc;

void Test_CreateWorm();

int mc::TestMarkov()
{
    sput_start_testing();
    sput_enter_suite("Test Updates:");
    sput_run_test(Test_CreateWorm);
    sput_finish_testing();
    return sput_get_return_value();
}

void Test_CreateWorm()
{
    para::ParaMC Para;
    Para.SetTest();
    weight::Weight Weight;
    Weight.SetTest(Para);
    diag::Diagram Diag;
    Diag.SetTest(Para.Lat, Para.RNG, Weight.G, Weight.W);
    Markov markov;
    markov.BuildNew(Para, Diag, Weight);

    int total = 0;
    for (int i = 0; i < 100; i++) {
        markov.CreateWorm();
        if (markov.Diag->Worm.Exist) {
            total += 1;
            //            markov.Diag->Save("diagram_template.config", "a");
        }
    }
    LOG_INFO("Update(Create Worm) are done!");
    sput_fail_unless(total == 100, "The accept ratio of CreateWorm = 1.0");

    Para.RNG.Reset(15);
    total = 0;
    for (int i = 0; i < 100; i++) {
        markov.CreateWorm();
        markov.DeleteWorm();
        if (!markov.Diag->Worm.Exist) {
            total += 1;
            //            markov.Diag->Save("diagram_template.config", "a");
        }
    }
    LOG_INFO("Update(Delete Worm) are done!");
    sput_fail_unless(total == 28, "The accept ratio of DeleteWorm = 0.25");
}
