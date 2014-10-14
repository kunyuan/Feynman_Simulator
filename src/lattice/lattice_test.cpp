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
    return sput_get_return_value();
}

void Test_Site()
{
    Site s;
    Vec<real> vec;
    vec[0]=5.5;
    vec[1]=7.5;
    
    s.Coordinate[0]=5;
    s.Coordinate[1]=7;
    s.SubLattice=1;
    sput_fail_unless(s.GetName()==235, "name of site(B, 5, 5) on 16*16");
    sput_fail_unless(s.GetVec()==vec, "vector of site(B, 5, 5) on 16*16");
}

