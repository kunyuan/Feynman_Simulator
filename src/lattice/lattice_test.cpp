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

int TestLattice()
{
    lattice.PlotLattice();
    
    sput_start_testing();
    sput_enter_suite("Test Definition of Class Site");
    sput_run_test(Test_Site);
    sput_enter_suite("Test Definition of Class Distance");
    sput_run_test(Test_Distance);
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
    sput_fail_unless(s.GetName()==235, "name of site(B, 5, 5) on 16*32");
    sput_fail_unless(s.GetVec()==vec, "vector of site(B, 5, 5) on 16*32");
    sput_fail_unless(s==GetSite(235), "Site of name 235");
}

void Test_Distance()
{
    Site s1(1), s2(0);
    
    s1.Coordinate[0] = 3;
    s1.Coordinate[1] = 1;
    
    s2.Coordinate[0] = 1;
    s2.Coordinate[1] = 31;
    
    Distance dis(s1, s2);   //s1---->s2, equal to dis=s2-s1;
    
    sput_fail_unless(dis.GetVec()==s2.GetVec()-s1.GetVec(), "distance=site1-site2, vecotr");
    sput_fail_unless(dis.SubLattice[0]==s1.SubLattice, "distance=site1-site2, sublattices");
    sput_fail_unless(dis.SubLattice[1]==s2.SubLattice, "distance=site1-site2, sublattices");
    sput_fail_unless(dis.Mirror().Dr[0]==2 && dis.Mirror().Dr[1]==2, "mirror(distance)");
}

