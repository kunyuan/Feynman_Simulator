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
using std::endl;

void Test_Site();
void Test_Distance();
void Test_Lattice();

Lattice lattice;

int TestLattice()
{
    lattice.PlotLattice();

    sput_start_testing();
    sput_enter_suite("Test Definition of Class Site");
    sput_run_test(Test_Site);
    sput_enter_suite("Test Definition of Class Distance");
    sput_run_test(Test_Distance);
    sput_enter_suite("Test Definition of Class Lattice");
    sput_run_test(Test_Lattice);
    sput_finish_testing();
    return sput_get_return_value();
}

void Test_Site()
{
    Site s;
    Vec<real> vec;
    if (D == 2) {
        vec[0] = 5.5;
        vec[1] = 7.5;

        s.Coordinate[0] = 5;
        s.Coordinate[1] = 7;
        s.Sublattice = 1;
        sput_fail_unless(s.GetName() == NSublattice * (L[0] * 7 + 5) + 1, "name of site(B, 5, 5)");
    }
}

void Test_Distance()
{
    Site s1(1), s2(0);

    if (D == 2) {
        s1.Coordinate[0] = 3;
        s1.Coordinate[1] = 1;

        s2.Coordinate[0] = 4;
        s2.Coordinate[1] = 3;

        Distance dis(s1, s2); //s1---->s2, equal to dis=s2-s1;

        sput_fail_unless(dis.Sublattice(IN) == s1.Sublattice, "distance=site1-site2, sublattices");
        sput_fail_unless(dis.Sublattice(OUT) == s2.Sublattice, "distance=site1-site2, sublattices");
        sput_fail_unless(dis.CoordiIndex() == 2 * L[0] + 1, "distance=site1-site2, sublattices");
        sput_fail_unless(dis.Mirror().Coordinate()[0] == 1 && dis.Mirror().Coordinate()[1] == 2, "mirror(distance)");
    }
}

void Test_Lattice()
{
    Site s(1);
    Vec<real> vec;
    Site s1(1), s2(0);

    if (D == 2) {
        vec[0] = 5.5;
        vec[1] = 7.5;

        s.Coordinate[0] = 5;
        s.Coordinate[1] = 7;

        s1.Coordinate[0] = 3;
        s1.Coordinate[1] = 1;

        s2.Coordinate[0] = 4;
        s2.Coordinate[1] = 3;

        Distance dis(s1, s2); //s1---->s2, equal to dis=s2-s1;
        sput_fail_unless(lattice.GetRealVec(s) == vec, "real vector of site(B, 5, 5)");
        sput_fail_unless(lattice.GetRealVec(dis) == lattice.GetRealVec(s2) - lattice.GetRealVec(s1), "real vector of distance");
    }
}
