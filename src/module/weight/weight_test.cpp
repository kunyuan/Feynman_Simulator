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
void TestWeightMeasuring();
void Sample(Sigma &, int num, Site in, Site out,
            spin SpinIn, spin SpinOut, real Beta);

int weight::TestWeight()
{
    sput_start_testing();
    sput_enter_suite("Test Weight...");

    //test diagram object weight, like sigma, G
    sput_run_test(TestWeightMeasuring);
    sput_finish_testing();
    return sput_get_return_value();
}

void TestWeightMeasuring()
{
    //some initialization
    Lattice lat(Vec<int>(1));
    real Beta = 1.0;
    weight::Sigma Sig(lat, Beta, 4);
    weight::Sigma Sig2(lat, Beta, 4);
    Site s1 = Site(0, 0);
    Site s2 = Site(0, 0);
    spin SpinIn = DOWN, SpinOut = DOWN;

    int Num = 1000000;
    //measure norm when set order=0
    Sample(Sig, Num, s1, s2, SpinIn, SpinOut, Beta);

    LOG_INFO("Order 1: " << Sig.WeightWithError(1) << endl
                         << "Order 2: " << Sig.WeightWithError(2) << endl
                         << "Order 3: " << Sig.WeightWithError(3));
    int order = Sig.OrderAcceptable(1, 500.0);
    LOG_INFO("Accepted Order=" << order);
    sput_fail_unless(order == 2, "Accepted order check.");
    Sig.UpdateWeight(2);
    LOG_INFO(Sig.Weight(s1, s2, 0.0, Beta / 2, SpinIn, SpinOut));
    //since order=1 and order=2 are accepted, so the weight will be
    //(P[1]+P[2])/P[0]*(average of each sample=0.5)=5/8

    //Weight class IO operation
    Sig.Save("test_weight.npz", "w");
    Sig2.Load("test_weight.npz");
    sput_fail_unless(Equal(Sig.Weight(s1, s2, 0.0, Beta / 2, SpinIn, SpinOut),
                           Sig2.Weight(s1, s2, 0.0, Beta / 2, SpinIn, SpinOut)),
                     "Weight class IO check.");
    system("rm test_weight.npz");
}

void Sample(Sigma &sigma, int num, Site in, Site out,
            spin SpinIn, spin SpinOut, real Beta)
{
    real P[4] = {16.0, 16.0, 4.0, 1.0};
    real PSum = P[0] + P[1] + P[2] + P[3];
    RandomFactory rng(100);
    for (int i = 0; i < num; i++) {
        if (i % 10 == 0)
            sigma.AddStatistics();
        real x = rng.urn() * PSum;
        if (x <= P[0])
            sigma.MeasureNorm();
        else if (x > P[0] && x <= P[0] + P[1])
            sigma.Measure(in, out, 0.0, rng.urn() * Beta,
                          SpinIn, SpinOut, 1, Complex(rng.urn(), rng.urn()));
        else if (x > P[0] + P[1] && x <= P[0] + P[1] + P[2])
            sigma.Measure(in, out, 0.0, rng.urn() * Beta,
                          SpinIn, SpinOut, 2, Complex(rng.urn(), rng.urn()));
        else
            sigma.Measure(in, out, 0.0, rng.urn() * Beta,
                          SpinIn, SpinOut, 3, Complex(rng.urn(), rng.urn()));
    }
}
