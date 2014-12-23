//
//  serialization_test.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 12/17/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//
#include <initializer_list>
#include <iostream>
#include <limits>
#include "utility/sput.h"
#include "utility/complex.h"
#include "dictionary.h"

using namespace std;

AnyObject ii = std::numeric_limits<int>::max();
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
    int IntMax = std::numeric_limits<int>::max();
    int IntMin = std::numeric_limits<int>::min();
    Port.Set("IntMax", IntMax);
    Port.Set("IntMin", IntMin);
    long long ago = std::numeric_limits<long long>::max();
    Port.Set("ago", ago);
    long long after = std::numeric_limits<long long>::min();
    Port.Set("after", after);
    unsigned long long biggest = std::numeric_limits<unsigned long long>::max();
    cout << biggest << endl;
    Port.Set("biggest", biggest);
    vector<int> v = { 1, 2, 3 };
    Port.Set("Vec", v);
    sput_fail_unless(Port.Get<int>("IntMax") == IntMax, "check integer type");
    sput_fail_unless(Port.Get<int>("IntMin") == IntMin, "check integer type");
    sput_fail_unless(Port.Get<long long>("ago") == ago,
                     "check long long integer type");
    sput_fail_unless(Port.Get<long long>("after") == after,
                     "check long long integer type");
    sput_fail_unless(Port.Get<unsigned long long>("biggest") == biggest,
                     "check unsigned long long integer type");
    auto vect = Port.Get<vector<int> >("Vec");
    sput_fail_unless(std::equal(v.begin(), v.end(),
                                vect.begin()),
                     "check vector<int> type");
    Complex ca = { 1.0, 2.0 };
    Complex cb = { 4.0, 2.1 };
    vector<Complex> vc = { ca, cb };
    Port.Set("cVec", vc);
    sput_fail_unless(Equal((Port.Get<vector<Complex> >("cVec"))[1], cb),
                     "check vector<Complex> type");
    Dictionary SubPort;
    Port.Name = "SubPort";
    SubPort.LoadFromString("{'b':11,'c':22}");
    Port.Set("dict", SubPort);
    sput_fail_unless(Port.Get<Dictionary>("dict").Get<int>("b") == 11,
                     "check dict type");
    Port.Save("test.txt", "w");
    Port.Clear();
    Port.Load("test.txt");
    sput_fail_unless(Equal((Port.Get<vector<Complex> >("cVec"))[1], cb),
                     "check vector<Complex> type");
    sput_fail_unless(Port.Get<Dictionary>("dict").Get<int>("b") == 11,
                     "check dict IO");
    //    Port.Print();
    system("rm test.txt");
    //    SubPort.Clear();
    //    SubPort = Port.Get<Dictionary>("dict");
    //    SubPort.Print();
}