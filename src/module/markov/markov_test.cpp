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
#include "module/weight/weight.h"
#include "module/parameter/parameter.h"
using namespace std;
using namespace mc;

void Test_Updates();

int mc::TestMarkov()
{
    sput_start_testing();
    sput_enter_suite("Test Updates:");
    sput_run_test(Test_Updates);
    sput_finish_testing();
    return sput_get_return_value();
}

void Test_Updates()
{
    para::ParaMC Para;
    Para.SetTest();
    weight::Weight Weight(true);
    Weight.SetTest(Para);
    diag::Diagram Diag;
    Diag.SetTest(Para.Lat, Para.RNG, Weight.G, Weight.W);
    Markov markov;
    markov.BuildNew(Para, Diag, Weight);

    //    Para.RNG.Reset(100);
    system("mkdir diagram");
    for (int i = 0; i < 100; i++) {
        //        if(i==4237)
        //            cout  << i << endl;
        markov.Hop(10000);
        sput_fail_unless(markov.Diag->CheckDiagram(), "Check for all the random steps");
        Diag.WriteDiagram2gv("diagram/" + ToString(i) + ".gv");
    }
    LOG_INFO("Updates Check are done!");
}
