//
//  markov_diag_counter.cpp
//  Feynman_Simulator
//
//  Created by yuan on 11/14/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include <stdio.h>
#include "markov.h"
#include "utility/sput.h"
#include "module/diagram/diagram.h"
#include "module/weight/weight.h"
#include "module/parameter/parameter.h"
using namespace std;
using namespace mc;

void Test_Counter();

int mc::TestDiagCounter()
{
    sput_start_testing();
    sput_enter_suite("Test Updates:");
    sput_run_test(Test_Counter);
    sput_finish_testing();
    return sput_get_return_value();
}

void Test_Counter()
{
    para::ParaMC Para;
    Para.SetTest();
    weight::Weight Weight(true);
    Weight.SetDiagCounter(Para);
    diag::Diagram Diag;
    Diag.SetTest(Para.Lat, Para.RNG, Weight.G, Weight.W);
    Markov markov;
    markov.BuildNew(Para, Diag, Weight);

//        Para.RNG.Reset(100);
    system("mkdir diagram");
    sput_fail_unless(Diag.CheckDiagram(), "Check diagram G,W,Ver and Weight");
    sput_fail_if(Equal(Diag.Weight,Complex(0.0, 0.0)), "Initialize diagram has no weight");
    int total=0;
    for (int i = 0; i < 1000; i++) {
        //        if(i==4237)
        //            cout  << i << endl;
        markov.Hop(100);
        
        if(Diag.Order==1)  total ++;
        sput_fail_unless(markov.Diag->CheckDiagram(), "Check for all the random steps");
        Diag.WriteDiagram2gv("diagram/" + ToString(i) + ".gv");
    }
    cout << "Number of Order1 diagram in 1000 samples:" << total <<endl;
    LOG_INFO("Updates Check are done!");
}
