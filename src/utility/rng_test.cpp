//
//  rng_test.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/5/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "rng.h"
#include "logger.h"
#include "timer.h"
#include "sput.h"
#include <string>
#include <sstream>
using namespace std;

void Test_RNG_IO()
{
    RandomFactory RNG;
    stringstream RngStr;
    RngStr << RNG;
    real a=RNG.urn();
    for(int i=0;i<100;i++)
        RNG.urn();
    RngStr >> RNG;
    sput_fail_unless(RNG.urn() == a, "import/export the RNG state");
}

void Test_RNG_Bound_And_Efficiency()
{
    const int bound = 5;
    int N = 9999999;
    bool flag = false;
    int Temp;
    double bin[bound] = {0};
    RandomFactory RNG;

    LOG_INFO("Real random number generator started..." << endl);
    timer T;
    T.start();

    for (int i = 0; i < N; i++) {
        RNG.urn();
    }
    T.stop();
    LOG_INFO("Time for " << N << " real numbers: " << T << endl);

    LOG_INFO("Int random number generator 0 started..." << endl);
    T.restart();
    for (int i = 0; i < N; i++) {
        Temp = RNG.irn(0, bound - 1);
        if (Temp < bound && Temp >= 0)
            bin[Temp]++;
        else
            flag = true;
    }
    T.stop();
    LOG_INFO("Time for " << N << " real numbers: " << T << endl);
    LOG_INFO("Distribution:" << endl);
    stringstream msg;
    for (int i = 0; i < bound; i++) {
        msg<<"Number=" << i << ", Prob=" << setprecision(4) << bin[i] / N << endl;
        bin[i] = 0.0;
    }
    LOG_INFO(msg.str());
    sput_fail_unless(flag == false, "RNG irn() exceeds the limit");

    flag = false;
    LOG_INFO("Int random number generator 1 started..." << endl);
    T.restart();
    for (int i = 0; i < N; i++) {
        Temp = RNG.irn1(0, bound - 1);
        if (Temp < bound && Temp >= 0)
            bin[Temp]++;
        else
            flag = true;
    }
    T.stop();
    LOG_INFO("Time for " << N << " real numbers: " << T << endl);
    LOG_INFO("Distribution:" << endl);
    msg.clear();
    for (int i = 0; i < bound; i++) {
        msg<<"Number=" << i << ", Prob=" << setprecision(4) << bin[i] / N << endl;
        bin[i] = 0.0;
    }
    LOG_INFO(msg.str());
    sput_fail_unless(flag == false, "RNG irn1() exceeds the limit");
}

int TestRNG()
{
    sput_start_testing();
    sput_enter_suite("Test Random Number Generator");
    sput_run_test(Test_RNG_IO);
    sput_run_test(Test_RNG_Bound_And_Efficiency);
    sput_finish_testing();
    return sput_get_return_value();
}
