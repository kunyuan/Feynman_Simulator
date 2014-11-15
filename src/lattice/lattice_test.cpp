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

using std::cout;
using std::endl;

void Test_Lattice();

int L[] = {16, 32};
Vec<int> size(L);

Lattice lattice(size, CHECKBOARD);

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
    if (lattice.Dimension == 2) {
        Vec<int> v0{1, 2};
        sput_fail_unless(v0[0] == 1 && v0[1] == 2, "Vector:test initialization");

        Vec<int> v, v1, v2;
        v[0] = 2;
        v[1] = 1;

        v2[0] = 8;
        v2[1] = 4;
        sput_fail_unless((v * 3 + v) == v2, "Vector:test the operators for vectors");

        Site s, s1, s2;
        Vec<real> vec;

        v[0] = 5;
        v[1] = 7;

        vec[0] = 5.5;
        vec[1] = 7.5;

        s.Coordinate = v;
        s.Sublattice = 1;

        sput_fail_unless(lattice.GetName(s) == NSublattice * (lattice.Size[1] * 5 + 7) + 1, "Site: name of site(B, 5, 7)");

        s1.Coordinate[0] = 3;
        s1.Coordinate[1] = 1;
        s1.Sublattice = 1;

        s2.Coordinate[0] = 4;
        s2.Coordinate[1] = 3;
        s2.Sublattice = 0;

        Distance dis = lattice.Dist(s1, s2); //s1---->s2, equal to dis=s2-s1;

        sput_fail_unless(lattice.GetSite(dis, IN) == s1, "Distance: distance=site1-site2, sublattices");
        sput_fail_unless(lattice.GetSite(dis, OUT) == s2, "Distance: distance=site1-site2, sublattices");
        sput_fail_unless(dis.CoordiIndex == 1 * lattice.Size[1] + 2, "Distance: distance=site1-site2, Coordinates");

        sput_fail_unless(lattice.GetRealVec(s) == vec, "Lattice: real vector for Site");
        sput_fail_unless(lattice.GetRealVec(dis) == (lattice.GetRealVec(s2) - lattice.GetRealVec(s1)), "Lattice: real vector of distance");
    }
    else if (lattice.Dimension == 3) {
        Vec<int> v, v1, v2;
        v[0] = 2;
        v[1] = 1;
        v[2] = 3;

        v2[0] = 8;
        v2[1] = 4;
        v2[2] = 12;
        sput_fail_unless((v * 3 + v) == v2, "Vector:test the operators for vectors");

        Site s, s1, s2;
        Vec<real> vec;

        v[0] = 5;
        v[1] = 7;
        v[2] = 2;

        vec[0] = 5.5;
        vec[1] = 7.5;
        vec[2] = 2.5;

        s.Coordinate = v;
        s.Sublattice = 1;

        sput_fail_unless(lattice.GetName(s) == NSublattice * (lattice.Size[0] * lattice.Size[1] * 2 + lattice.Size[0] * 7 + 5) + 1, "Site: name of site(B, 5, 7, 2)");

        s1.Coordinate[0] = 3;
        s1.Coordinate[1] = 1;
        s1.Coordinate[2] = 4;
        s1.Sublattice = 1;

        s2.Coordinate[0] = 4;
        s2.Coordinate[1] = 3;
        s2.Coordinate[2] = 7;
        s2.Sublattice = 0;

        Distance dis = lattice.Dist(s1, s2); //s1---->s2, equal to dis=s2-s1;

        sput_fail_unless(lattice.GetSite(dis, IN) == s1, "Distance: distance=site1-site2, sublattices");
        sput_fail_unless(lattice.GetSite(dis, OUT) == s2, "Distance: distance=site1-site2, sublattices");
        sput_fail_unless(dis.CoordiIndex == 3 * lattice.Size[0] * lattice.Size[1] + 2 * lattice.Size[0] + 1, "Distance: distance=site1-site2, Coordinates");

        sput_fail_unless(lattice.GetRealVec(s) == vec, "Lattice: real vector for Site");
        sput_fail_unless(lattice.GetRealVec(dis) == (lattice.GetRealVec(s2) - lattice.GetRealVec(s1)), "Lattice: real vector of distance");
    }
}
