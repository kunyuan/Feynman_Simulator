//
//  observable_test.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/13/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "weight.h"
#include "component.h"
#include "utility/sput.h"
#include "utility/rng.h"
#include "utility/dictionary.h"

using namespace std;
using namespace weight;

void TestDiagramObject();
void TestWeightMeasuring();
void WeightMeasuring(real Beta, int Num);
void Sample(Sigma&, int num, Site in, Site out,
            spin SpinIn, spin SpinOut, real Beta, int order, real* P);

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
    real Beta = 5.0;
    int num = 100000;
    WeightMeasuring(Beta, num);
    WeightMeasuring(Beta, num * 10);
    WeightMeasuring(Beta / 10, num);
}

void WeightMeasuring(real Beta, int Num)
{
    //some initialization
    Lattice lat(Vec<int>(4));
    weight::Sigma Sig(lat, Beta, 32, 4);
    weight::Sigma Sig2(lat, Beta, 32, 4);
    Site s1 = Site(0, 0);
    Site s2 = Site(0, 0);
    spin SpinIn = DOWN, SpinOut = DOWN;
    real tau = Beta / 2;
    int order = 4;

    //provide order+1 real numbers, order=0 corresponds to the Norm term
    real P[] = { 16.0, 16.0, 4.0, 1.0, 1.0 };
    //measure norm when set order=0
    Sample(Sig, Num, s1, s2, SpinIn, SpinOut, Beta, order + 1, P);

    LOG_INFO("Order 1: " << Sig.Estimator.RelativeError(1) << endl
                         << "Order 2: " << Sig.Estimator.RelativeError(2) << endl
                         << "Order 3: " << Sig.Estimator.RelativeError(3));
    LOG_INFO("Accepted Order with Threshold 0.04=" << Sig.Estimator.OrderAcceptable(1, 0.04));

    Sig.UpdateWeight(2); //accept up to order 2
    Complex weight = Sig.Weight(s1, s2, 0.0, tau, SpinIn, SpinOut);
    real w = (P[1] + P[2]) / P[0] * 0.5 / Beta; //0.5 is the average <rng.urn()>
    Complex realweight = Complex(w, w);
    //2 standard deviation is allowd in the estimate
    Complex relative_error = 2.0 * (Sig.Estimator.RelativeError(1) + Sig.Estimator.RelativeError(2));
    Complex error(realweight.Re * relative_error.Re,
                  realweight.Im * relative_error.Im);

    LOG_INFO(weight);
    sput_fail_unless(Equal(weight.Re, realweight.Re, error.Re),
                     "Check real part of Weight and it's error");
    sput_fail_unless(Equal(weight.Im, realweight.Im, error.Im),
                     "Check imag part of Weight and it's error");

    //Weight class IO operation

    Dictionary dict = Sig.ToDict();
    Sig2.FromDict(dict);
    sput_fail_unless(Equal(Sig.Weight(s1, s2, 0.0, Beta / 2, SpinIn, SpinOut),
                           Sig2.Weight(s1, s2, 0.0, Beta / 2, SpinIn, SpinOut)),
                     "Weight class IO check.");
}

void Sample(Sigma& sigma, int num, Site in, Site out,
            spin SpinIn, spin SpinOut, real Beta, int order, real* P)
{
    RandomFactory rng(100);
    real P_lower[order + 1];
    real P_upper[order + 1];
    P_lower[0] = 0.0;
    P_upper[0] = P[0];
    for (int i = 0; i < order; i++) {
        P_lower[i + 1] = P_upper[i];
        P_upper[i + 1] = P_lower[i + 1] + P[i + 1];
    }
    real PSum = P_upper[order];

    for (int i = 0; i < num; i++) {
        if (i % 100 == 0)
            sigma.Estimator.AddStatistics();
        real x = rng.urn() * PSum;
        for (int o = 0; o <= order; o++) {
            if (x > P_lower[o] && x <= P_upper[o]) {
                if (o == 0)
                    sigma.Estimator.MeasureNorm();
                else
                    sigma.Measure(in, out, 0.0, rng.urn() * Beta,
                                  SpinIn, SpinOut, o, Complex(rng.urn(), rng.urn()));
                break;
            }
        }
    }
}
