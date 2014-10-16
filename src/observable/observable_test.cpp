//
//  observable_test.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/13/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "observable.h"
#include "cnpy.h"
#include "sput.h"
using namespace std;

void TestObservableIO();

int TestObservable()
{
    sput_start_testing();
    sput_enter_suite("Test Observable...");
    sput_run_test(TestObservableIO);
    sput_finish_testing();
    return sput_get_return_value();
}

void TestObservableIO()
{
    EstimateKeeper<Complex> quan("TestQuan");
    Complex a[10];
    for(int i=0;i<10;i++)
    {
        a[i]=Complex(i+1,10-i);
        quan.AddStatistics(a[i]);
    }
    Estimate<Complex> ExpectedResult(Complex(5.0, 5.0),Complex(1.5, 0.9));
    //!!!This two value only works if you set _norm=1.0 and ThrowRatio=1.0/3
    sput_fail_unless(Equal(quan.Estimate().Mean, ExpectedResult.Mean),
                     "check the Mean value.");
    sput_fail_unless(Equal(quan.Estimate().Error, ExpectedResult.Error),
                     "check the Error value.");
    
    //EstimateKeeper IO operation
    quan.WriteToFile("TestObservable", "w");
    cnpy::npz_t test=cnpy::npz_load("TestObservable.npz");
    EstimateKeeper<Complex> quan_another("TestQuan");
    quan_another.ReadFromFile(test);
    sput_fail_unless(Equal(quan_another.Estimate().Mean, ExpectedResult.Mean),
                     "check the Mean value.");
    sput_fail_unless(Equal(quan_another.Estimate().Error, ExpectedResult.Error),
                     "check the Error value.");
}