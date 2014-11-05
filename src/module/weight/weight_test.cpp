//
//  observable_test.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/13/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "weight.h"
#include "utility/sput.h"
#include "utility/rng.h"

using namespace std;
using namespace weight;

void TestDiagramObject();

int weight::TestWeight()
{
    sput_start_testing();
    sput_enter_suite("Test Weight...");

    //test diagram object weight, like sigma, G
    sput_run_test(TestDiagramObject);
    sput_finish_testing();
    return sput_get_return_value();
}

void TestDiagramObject()
{
    RandomFactory rng(100);
    Lattice lat(Vec<int>(1));
    real Beta = 5.0;
    weight::Sigma Sig(lat, Beta, 4);
    Site s1, s2;
    for (int i = 0; i < 1000000; i++) {
        Sig.Measure(s1, s2, 0.0, rng.urn() * Beta, DOWN, DOWN, 1, Complex(rng.urn(), rng.urn()));
        if (i % 4 == 0)
            Sig.Measure(s1, s2, 0.0, rng.urn() * Beta, DOWN, DOWN, 2, Complex(rng.urn(), rng.urn()));
        if (i % 9 == 0) {
            Sig.Measure(s1, s2, 0.0, rng.urn() * Beta, DOWN, DOWN, 3, Complex(rng.urn(), rng.urn()));
            Sig.AddStatistics();
        }
    }
    //Do random walk so that the error ratio between Order 3, Order 2 and Order 1 is 3:2:1 roughly
    /*
      After the random walk, the first order of Sig has the weight mean/Beta/(1+1/4+1/9),
      where mean=0.5 is the mean of samples,
      the strange factor in denominator comes the fact that Sig._Norm
      will increase for all order 1, 2 and 3 measuring. 
    */

    LOG_INFO("Order 1: " << Sig.WeightWithError(1) << endl << "Order 2: " << Sig.WeightWithError(2) << endl << "Order 3: " << Sig.WeightWithError(3));
    int order = Sig.OrderAcceptable(1, 500.0);
    LOG_INFO("Accepted Order=" << order);
    sput_fail_unless(order == 2, "Accepted order check.");
    Sig.UpdateWeight(2);
    LOG_INFO(Sig.Weight(s1, s2, 0.0, Beta / 2, DOWN, DOWN));

    //Weight class IO operation
    Sig.Save("test_weight.npz", "w");
    weight::Sigma Sig2(lat, Beta, 4);
    Sig2.Load("test_weight.npz");
    sput_fail_unless(Equal(Sig.Weight(s1, s2, 0.0, Beta / 2, DOWN, DOWN), Sig.Weight(s1, s2, 0.0, Beta / 2, DOWN, DOWN)), "Weight class IO check.");
    system("rm test_weight.npz");
}