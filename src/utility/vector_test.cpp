//
//  vector_test.cpp
//  Feynman_Simulator
//
//  Created by yuan on 10/14/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include <iostream>
#include "vector.h"
#include "sput.h"

using std::cout;
using std::endl;

void Test_Vector();

int TestVector()
{
    sput_start_testing();
    sput_enter_suite("Test Definition of Vector");
    sput_run_test(Test_Vector);
    sput_finish_testing();
    return sput_get_return_value();
}

void Test_Vector()
{
    Vec<int> v, v2;
    v[0]=2;
    v[1]=1;
    v2[0]=8;
    v2[1]=4;
    cout << v2 << ' ' << v << v.PrettyString() <<endl;
    sput_fail_unless((v*3+v)==v2, "test the operators for vectors");
}

