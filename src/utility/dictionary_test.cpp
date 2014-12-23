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
using namespace Python;

void Test_Ref();
void Test_Cast();
void Test_Dict();
int TestDictionary()
{
    sput_start_testing();
    sput_enter_suite("Test Python wrapper layer");
    sput_run_test(Test_Ref);
    sput_run_test(Test_Cast);
    sput_run_test(Test_Dict);
    sput_finish_testing();
    return sput_get_return_value();
}

void Test_Ref()
{
    PyObject* iptr = PyInt_FromLong(1);
    auto i_ref = iptr->ob_refcnt;
    AnyObject i1 = iptr;
    sput_fail_unless(i1.RefCount() == i_ref, "create AnyObject from PyObject* who owns the reference");
    AnyObject i2 = Object(iptr, NoRef);
    sput_fail_unless(i2.RefCount() == i_ref + 1, "create AnyObject from PyObject* who does not own the reference");
    AnyObject i3 = i2;
    sput_fail_unless(i3.RefCount() == i_ref + 2, "create AnyObject with copy constructor");
    i2.Destroy();

    AnyObject j = 2;
    auto j_ref = j.RefCount();
    i3 = j;
    //all new reference to iptr has been destroyed except i1
    sput_fail_unless(iptr->ob_refcnt == i_ref, "the old AnyObject has to call Py_DEREF");
    sput_fail_unless(i3.RefCount() == j_ref + 1, "the new AnyObject own an new reference");
    i3.Destroy();

    PyObject* jptr = j.Get(NewRef);
    sput_fail_unless(jptr->ob_refcnt == j_ref + 1, "get PyObject* with new reference");
    Py_XDECREF(jptr);
    PyObject* jjptr = j.Get(NoRef);
    sput_fail_unless(jjptr->ob_refcnt == j_ref, "get PyObject* witout reference");
}

void Test_Cast()
{
    AnyObject value;
    int IntMax = std::numeric_limits<int>::max();
    value = IntMax;
    sput_fail_unless(IntMax == value.As<int>(), "Int Max");
    int IntMin = std::numeric_limits<int>::min();
    value = IntMin;
    sput_fail_unless(IntMin == value.As<int>(), "Int Min");

    long long lliMax = std::numeric_limits<long long>::max();
    value = lliMax;
    sput_fail_unless(lliMax == value.As<long long>(), "Long Long Int Max");
    long long lliMin = std::numeric_limits<long long>::min();
    value = lliMin;
    sput_fail_unless(lliMin == value.As<long long>(), "Long Long Int Min");
    unsigned long long biggest = std::numeric_limits<unsigned long long>::max();
    value = biggest;
    sput_fail_unless(biggest == value.As<unsigned long long>(),
                     "Unsigned Long Long Int Max");
}

void Test_Dict()
{
    Dictionary Port;
    unsigned long long biggest = std::numeric_limits<unsigned long long>::max();
    Port.Set("biggest", biggest);
    vector<int> v = { 1, 2, 3 };
    Port.Set("Vec", v);
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
    SubPort.Clear();
    SubPort = Port.Get<Dictionary>("dict");
    SubPort.Print();

    sput_fail_unless(Port.RefCount() == 1, "Port ref check");
    Port.Destroy();
    sput_fail_unless(SubPort.RefCount() == 1, "SubPort ref check");
    SubPort.Destroy();
}