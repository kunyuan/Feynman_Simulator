//
//  lattice_test.cpp
//  Fermion_Simulator
//
//  Created by yuan on 10/10/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include <iostream>
#include "lattice.h"
#include "sput.h"

using std::cout;

void Test_Site();
void Test_Distance();

int TestLattice()
{
    sput_start_testing();
    sput_enter_suite("Test Definition of Class Site");
    sput_run_test(Test_Site);
    sput_finish_testing();
    a();
    return sput_get_return_value();
}

void Test_Site()
{
//    Site s;
}

