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
    s.Sublattice=1;
    sput_fail_unless(s.GetName()==NSublattice*(L[0]*7+5)+1, "name of site(B, 5, 5)");
    sput_fail_unless(s.GetVec()==vec, "vector of site(B, 5, 7) ");
}

void Test_Distance()
{
    Site s1(1), s2(0);
    
    s1.Coordinate[0] = 3;
    s1.Coordinate[1] = 1;
    
    s2.Coordinate[0] = 4;
    s2.Coordinate[1] = 3;
    
    Distance dis(s1, s2);   //s1---->s2, equal to dis=s2-s1;
    
    sput_fail_unless(dis.GetVec()==s2.GetVec()-s1.GetVec(), "distance=site1-site2, vecotr");
    sput_fail_unless(dis.GetSublattice(IN)==s1.Sublattice, "distance=site1-site2, sublattices");
    sput_fail_unless(dis.GetSublattice(OUT)==s2.Sublattice, "distance=site1-site2, sublattices");
    sput_fail_unless(dis.Coordinate()==2*L[0]+1, "distance=site1-site2, sublattices");
    sput_fail_unless(dis.Mirror().dCoordinate[0]==1 && dis.Mirror().dCoordinate[1]==2, "mirror(distance)");
}

