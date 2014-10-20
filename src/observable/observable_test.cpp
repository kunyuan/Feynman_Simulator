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
    QuanVector2.ReadState("TestObservable");
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
    RandomFactory rng;
    Lattice lat;
    real Beta = 5.0;
    Weight::Sigma Sig(lat, Beta);
    Distance d(0, 0);
    for (int i = 0; i < 1000000; i++) {
        Sig.Measure(Complex(rng.nrn(2.0, 1.0), rng.nrn(2.0, 1.0)), d, Beta * rng.urn(), UP, UP);
    }
    Sig.SaveState("test_weight");
}