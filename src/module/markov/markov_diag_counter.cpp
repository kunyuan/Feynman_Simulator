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
    Diag.SetTest(Para.Lat, *Weight.G, *Weight.W);
    Markov markov;
    markov.BuildNew(Para, Diag, Weight);

    system("rm -rf diagram");
    system("mkdir diagram");
    sput_fail_unless(Diag.CheckDiagram(), "Check diagram G,W,Ver and Weight");
    sput_fail_if(Equal(Diag.Weight, Complex(0.0, 0.0)), "Initialize diagram has nonzero weight");

    int total[MAX_ORDER] = { 0 };

    for (int i = 0; i < 5000; i++) {
        //        if(i==0)
        //            cout << i<<endl;
        markov.Hop(40000);

        sput_fail_unless(markov.Diag->CheckDiagram(), "Check for all the random steps");
        if (!markov.Diag->Worm.Exist) {
            total[Diag.Order]++;
            Diag.WriteDiagram2gv("diagram/" + ToString(Para.Counter) + ".gv");
        }
    }
    cout << "Number of different Order diagrams in 5000 samples: " << real(total[2]) / real(total[1]) << " " <<  4.0*real(total[3]) / real(total[1]) << endl;
    LOG_INFO("Updates Check are done!");
}
