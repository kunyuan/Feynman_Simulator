//
//  observable_test.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/13/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "weight.h"
#include "cnpy.h"
#include "sput.h"
#include "utility.h"
#include "convention.h"
#include "rng.h"

using namespace std;

void TestObservableComplex();
void TestObservableReal();
void TestDiagramObject();

int TestObservable()
{
    sput_start_testing();
    sput_enter_suite("Test Observable...");

    //complex and real Estimator have to be tested separatelly,
    //since Estimator._update() is different for those two types
    sput_run_test(TestObservableComplex);
    sput_run_test(TestObservableReal);

    //test diagram object weight, like sigma, G
    sput_run_test(TestDiagramObject);
    sput_finish_testing();
    return sput_get_return_value();
}

void TestObservableComplex()
{
    Estimator<Complex> quan1("1");
    Estimator<Complex> quan2("2");
    Complex a[10];
    for (int i = 0; i < 10; i++) {
        a[i] = Complex(i + 1, 10 - i);
        quan1.Measure(a[i]);
        quan1.AddStatistics();
        quan2.Measure(-a[i]);
        quan2.AddStatistics();
    }
    Estimate<Complex> ExpectedResult(Complex(5.0, 5.0), Complex(1.5, 0.9));
    //!!!This two value only works if you set _norm=1.0 and ThrowRatio=1.0/3
    sput_fail_unless(Equal(quan1.Estimate().Mean, ExpectedResult.Mean),
                     "check the Mean value.");
    sput_fail_unless(Equal(quan1.Estimate().Error, ExpectedResult.Error),
                     "check the Error value.");

    //Estimator IO operation

    EstimatorBundle<Complex> QuanVector;
    QuanVector.AddEstimator(quan1);
    QuanVector.AddEstimator(quan2);
    QuanVector.SaveState("TestObservable", "w");
    EstimatorBundle<Complex> QuanVector2;
    QuanVector2.AddEstimator("1");
    QuanVector2.AddEstimator("2");
    QuanVector2.LoadState("TestObservable");
    sput_fail_unless(Equal(QuanVector2[1].Estimate().Mean, -ExpectedResult.Mean),
                     "EstimatorVector:check the Mean value.");
    sput_fail_unless(Equal(QuanVector2[1].Estimate().Error, ExpectedResult.Error),
                     "EstimatorVector:check the Error value.");
    sput_fail_unless(Equal(QuanVector2["2"].Estimate().Mean, -ExpectedResult.Mean),
                     "EstimatorVector:check the Mean value.");
    sput_fail_unless(Equal(QuanVector2["2"].Estimate().Error, ExpectedResult.Error),
                     "EstimatorVector:check the Error value.");
}

void TestObservableReal()
{
    Estimator<real> quan1("1");
    Estimator<real> quan2("2");
    real a[10];
    for (int i = 0; i < 10; i++) {
        a[i] = i + 1;
        quan1.Measure(a[i]);
        quan1.AddStatistics();
    }
    Estimate<Complex> ExpectedResult(5.0, 1.5);
    //!!!This two value only works if you set _norm=1.0 and ThrowRatio=1.0/3
    sput_fail_unless(Equal(quan1.Estimate().Mean, ExpectedResult.Mean),
                     "check the Mean value.");
    sput_fail_unless(Equal(quan1.Estimate().Error, ExpectedResult.Error),
                     "check the Error value.");
}

void TestDiagramObject()
{
    RandomFactory rng(100);
    Lattice lat;
    real Beta = 5.0;
    Weight::Sigma Sig(lat, Beta, 4);
    Distance d(0, 0);
    for (int i = 0; i < 1000000; i++) {
        Sig.Measure(d, rng.urn() * Beta, DOWN, DOWN, 1, Complex(rng.urn(), rng.urn()));
        if (i % 4 == 0)
            Sig.Measure(d, rng.urn() * Beta, DOWN, DOWN, 2, Complex(rng.urn(), rng.urn()));
        if (i % 9 == 0) {
            Sig.Measure(d, rng.urn() * Beta, DOWN, DOWN, 3, Complex(rng.urn(), rng.urn()));
            Sig.AddStatistics();
        }
    }
    //Do random walk so that the error ratio between Order 3, Order 2 and Order 1 is 3:2:1 roughly

    LOG_INFO("Order 1: " << Sig.WeightWithError(1) << endl << "Order 2: " << Sig.WeightWithError(2) << endl << "Order 3: " << Sig.WeightWithError(3));
    int order = Sig.OrderAcceptable(1, 500.0);
    LOG_INFO("Accepted Order=" << order);
    sput_fail_unless(order == 2, "Accepted order check.");
    Sig.UpdateWeight(2);
    LOG_INFO(Sig.Weight(d, Beta / 2, DOWN, DOWN));

    //Weight class IO operation
    Sig.SaveState("test_weight", "w");
    Weight::Sigma Sig2(lat, Beta, 4);
    Sig2.LoadState("test_weight");
    sput_fail_unless(Equal(Sig.Weight(d, Beta / 2, DOWN, DOWN), Sig.Weight(d, Beta / 2, DOWN, DOWN)), "Weight class IO check.");
}