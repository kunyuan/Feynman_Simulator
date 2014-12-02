//
//  dyson_test.cpp
//  Feynman_Simulator
//
//  Created by yuan on 11/19/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include <stdio.h>
#include "dyson.h"
#include "module/weight/weight.h"
#include "module/weight/weight_inherit.h"
#include "module/parameter/parameter.h"
#include "environment/environment.h"
#include "utility/sput.h"
#include "utility/utility.h"
#include "utility/rng.h"

using namespace std;
using namespace weight;
using namespace dyson;

void TestMultiply();
void TestInverse();
void TestG();
void TestW();

int dyson::TestDyson()
{
    sput_start_testing();
    sput_enter_suite("Test Dyson...");
    
    sput_run_test(TestMultiply);
    sput_run_test(TestInverse);
    sput_run_test(TestG);
    sput_run_test(TestW);

    sput_finish_testing();
    return sput_get_return_value();
}

void TestMultiply()
{
    Complex mat1[40];
    Complex mat2[40];
    for(int i=0; i<10; i++)
    {
        mat1[0*10+i] = Complex(1.0, 0.0);
        mat1[1*10+i] = Complex(2.0, 0.0);
        mat1[2*10+i] = Complex(2.0, 0.0);
        mat1[3*10+i] = Complex(3.0, 0.0);
        
        mat2[0*10+i] = Complex(3.0, 0.0);
        mat2[1*10+i] = Complex(2.0, 0.0);
        mat2[2*10+i] = Complex(2.0, 0.0);
        mat2[3*10+i] = Complex(1.0, 0.0);
    }
    Complex mat3[40];
    AssignFromTo(mat1, mat3, 40);
    MatrixMultiplySUB(mat3, mat2, 10);
    
    sput_fail_unless(Equal(mat3[0], Complex(7.0, 0.0)), "Check: matrix sub multiply");
    sput_fail_unless(Equal(mat3[5], Complex(7.0, 0.0)), "Check: matrix sub multiply");
    sput_fail_unless(Equal(mat3[10], Complex(4.0, 0.0)), "Check: matrix sub multiply");
    
    Complex mat4[4]={Complex(3.0, 0.0),
                    Complex(2.0, 0.0),
                    Complex(2.0, 0.0),
                    Complex(1.0, 0.0)};
    AssignFromTo(mat1, mat3, 40);
    MatrixMultiplySUB(mat3, mat4, 10, 1);
    
    sput_fail_unless(Equal(mat3[0], Complex(7.0, 0.0)), "Check: matrix sub multiply");
    sput_fail_unless(Equal(mat3[5], Complex(7.0, 0.0)), "Check: matrix sub multiply");
    sput_fail_unless(Equal(mat3[10], Complex(4.0, 0.0)), "Check: matrix sub multiply");
    sput_fail_unless(Equal(mat3[22], Complex(12.0, 0.0)), "Check: matrix sub multiply");
    
    Complex mat5[48]={0.0};
    Complex mat6[48]={0.0};
    for(int i=0; i<3; i++)
    {
        mat5[0*3+i] = Complex(1.0, 0.0);
        mat5[1*3+i] = Complex(-1.0, 0.0);
        mat5[4*3+i] = Complex(-1.0, 0.0);
        mat5[5*3+i] = Complex(1.0, 0.0);
        mat5[11*3+i] = Complex(2.0, 0.0);
        mat5[14*3+i] = Complex(2.0, 0.0);
        
        mat6[0*3+i] = Complex(1.0, 0.0);
        mat6[1*3+i] = Complex(1.0, 0.0);
        mat6[4*3+i] = Complex(1.0, 0.0);
        mat6[5*3+i] = Complex(1.0, 0.0);
        mat6[11*3+i] = Complex(2.0, 0.0);
        mat6[14*3+i] = Complex(2.0, 0.0);
    }
    
    Complex mat7[48];
    AssignFromTo(mat5, mat7, 48);
    MatrixMultiplySP(mat7, mat6, 3);
    
    sput_fail_unless(Equal(mat7[0], Complex(0.0, 0.0)), "Check: matrix sp multiply");
    sput_fail_unless(Equal(mat7[1*3], Complex(0.0, 0.0)), "Check: matrix sp multiply");
    sput_fail_unless(Equal(mat7[4*3], Complex(0.0, 0.0)), "Check: matrix sp multiply");
    sput_fail_unless(Equal(mat7[5*3], Complex(0.0, 0.0)), "Check: matrix sp multiply");
    sput_fail_unless(Equal(mat7[11*3], Complex(0.0, 0.0)), "Check: matrix sp multiply");
    sput_fail_unless(Equal(mat7[14*3], Complex(0.0, 0.0)), "Check: matrix sp multiply");
    sput_fail_unless(Equal(mat7[10*3], Complex(4.0, 0.0)), "Check: matrix sp multiply");
    sput_fail_unless(Equal(mat7[15*3], Complex(4.0, 0.0)), "Check: matrix sp multiply");
    sput_fail_unless(Equal(mat7[10*3+2], Complex(4.0, 0.0)), "Check: matrix sp multiply");
    sput_fail_unless(Equal(mat7[15*3+1], Complex(4.0, 0.0)), "Check: matrix sp multiply");
}

void TestInverse()
{
    Complex mat1[40];
    for(int i=0; i<10; i++)
    {
        mat1[0*10+i] = Complex(1.0, 0.0);
        mat1[1*10+i] = Complex(2.0, 0.0);
        mat1[2*10+i] = Complex(2.0, 0.0);
        mat1[3*10+i] = Complex(3.0, 0.0);
        
    }
    Complex mat2[40];
    AssignFromTo(mat1, mat2, 40);
    MatrixInverseSUB(mat2, 10);
    MatrixMultiplySUB(mat2, mat1, 10);
    
    sput_fail_unless(Equal(mat2[0], Complex(1.0, 0.0)), "Check: matrix sub multiply");
    sput_fail_unless(Equal(mat2[5], Complex(1.0, 0.0)), "Check: matrix sub multiply");
    sput_fail_unless(Equal(mat2[10], Complex(0.0, 0.0)), "Check: matrix sub multiply");
    sput_fail_unless(Equal(mat2[20], Complex(0.0, 0.0)), "Check: matrix sub multiply");
    sput_fail_unless(Equal(mat2[30], Complex(1.0, 0.0)), "Check: matrix sub multiply");
    
    Complex mat3[64]={0.0};
    for(int i=0; i<4; i++)
    {
        mat3[0*4+i] = Complex(0.0, 0.0);
        mat3[1*4+i] = Complex(1.0, 0.0);
        mat3[4*4+i] = Complex(1.0, 0.0);
        mat3[5*4+i] = Complex(0.0, 0.0);
        
        mat3[10*4+i] = Complex(0.0, 0.0);
        mat3[11*4+i] = Complex(2.0, 0.0);
        mat3[14*4+i] = Complex(2.0, 0.0);
        mat3[15*4+i] = Complex(0.0, 0.0);
        
    }
    Complex mat4[64];
    AssignFromTo(mat3, mat4, 64);
    MatrixInverseSP(mat4, 4);
    MatrixMultiplySP(mat4, mat3, 4);
    
    sput_fail_unless(Equal(mat4[0],      Complex(1.0, 0.0)), "Check: matrix multiply");
    sput_fail_unless(Equal(mat4[1*4+1],  Complex(0.0, 0.0)), "Check: matrix multiply");
    sput_fail_unless(Equal(mat4[4*4+2],  Complex(0.0, 0.0)), "Check: matrix multiply");
    sput_fail_unless(Equal(mat4[5*4+3],  Complex(1.0, 0.0)), "Check: matrix multiply");
    sput_fail_unless(Equal(mat4[10*4+2], Complex(1.0, 0.0)), "Check: matrix multiply");
    sput_fail_unless(Equal(mat4[11*4+2], Complex(0.0, 0.0)), "Check: matrix multiply");
    sput_fail_unless(Equal(mat4[15*4+1], Complex(1.0, 0.0)), "Check: matrix multiply");
}

void TestG()
{
    para::ParaDyson Para;
    Para.SetTest();
    weight::Weight Weight(true);
    Weight.BuildNew(GW|SigmaPolar,Para);
}

void TestW()
{
    
}
