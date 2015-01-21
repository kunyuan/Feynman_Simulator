//
//  lattice_test.cpp
//  Fermion_Simulator
//
//  Created by yuan on 10/10/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "lattice.h"
#include <iostream>
#include "../utility/sput.h"
#include "../utility/convention.h"

using std::cout;
using std::endl;

void Test_Lattice();

int L[] = { 16, 32 };
Vec<int> size(L);

Lattice lattice(size, 2);

int TestLattice()
{
    sput_start_testing();
    sput_enter_suite("Test Definition of Class Lattice");
    sput_run_test(Test_Lattice);
    sput_finish_testing();
    return sput_get_return_value();
}

void Test_Lattice()
{
    Vec<int> v0{ 1, 2 };
    sput_fail_unless(v0[0] == 1 && v0[1] == 2, "Vector:test initialization");

    Vec<int> v, v1, v2;
    v = { 2, 1 };
    v2 = { 8, 4 };
    sput_fail_unless((v * 3 + v) == v2, "Vector:test the operators for vectors");
    Site s1, s2;
    s1.Coordinate = { 3, 1 };
    s1.Sublattice = 1;
    s2.Coordinate = { 4, 3 };
    s2.Sublattice = 0;
}
