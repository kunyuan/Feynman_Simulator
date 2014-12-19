//
//  serialization_test.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 12/17/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//
#include <initializer_list>
#include <iostream>
#include "utility/sput.h"
#include "utility/complex.h"
#include "dictionary.h"

using namespace std;

void Test_Dict();
int TestDictionary()
{
    sput_start_testing();
    sput_enter_suite("Test Definition of Class Lattice");
    sput_run_test(Test_Dict);
    sput_finish_testing();
    return sput_get_return_value();
}

void Test_Dict()
{
    Dictionary Port;
    Port.Set("a", 1);
    vector<int> v = { 1, 2, 3 };
    Port.Set("Vec", v);
    sput_fail_unless(Port.Get<int>("a") == 1, "check integer type");
    sput_fail_unless(std::equal(v.begin(), v.end(),
                                Port.Get<vector<int> >("Vec").begin()),
                     "check vector<int> type");
    Complex a = { 1.0, 2.0 };
    Complex b = { 4.0, 2.1 };
    vector<Complex> vc = { a, b };
    Port.Set("cVec", vc);
    sput_fail_unless(Equal((Port.Get<vector<Complex> >("cVec"))[1], b),
                     "check vector<Complex> type");
    Dictionary SubPort;
    SubPort.LoadByEval("{'b':11,'c':22}");
    Port.Set("dict", SubPort);
    sput_fail_unless(Port.Get<Dictionary>("dict").Get<int>("b") == 11,
                     "check dict type");
}