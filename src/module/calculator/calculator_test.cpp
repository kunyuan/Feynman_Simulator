//
//  dyson_test.cpp
//  Feynman_Simulator
//
//  Created by yuan on 11/19/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "calculator.h"
#include "module/weight/weight.h"
#include "module/weight/weight_component.h"
#include "module/parameter/parameter.h"
#include "environment/environment.h"
#include "utility/sput.h"
#include "utility/utility.h"
#include "utility/rng.h"

using namespace std;
using namespace weight;
using namespace calc;

void TestMultiply();
void TestInverse();
void TestG();
void TestW();

int calc::TestCalculator()
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
    for (int i = 0; i < 10; i++) {
        mat1[0 * 10 + i] = Complex(1.0, 0.0);
        mat1[1 * 10 + i] = Complex(2.0, 0.0);
        mat1[2 * 10 + i] = Complex(2.0, 0.0);
        mat1[3 * 10 + i] = Complex(3.0, 0.0);

        mat2[0 * 10 + i] = Complex(3.0, 0.0);
        mat2[1 * 10 + i] = Complex(2.0, 0.0);
        mat2[2 * 10 + i] = Complex(2.0, 0.0);
        mat2[3 * 10 + i] = Complex(1.0, 0.0);
    }
    Complex mat3[40];
    AssignFromTo(mat1, mat3, 40);
    calc::MatrixMultiply(mat3, mat2, 10);

    sput_fail_unless(Equal(mat3[0], Complex(7.0, 0.0)), "Check: matrix multiply");
    sput_fail_unless(Equal(mat3[5], Complex(7.0, 0.0)), "Check: matrix multiply");
    sput_fail_unless(Equal(mat3[10], Complex(4.0, 0.0)), "Check: matrix multiply");

    Complex mat4[4] = {Complex(3.0, 0.0),
                       Complex(2.0, 0.0),
                       Complex(2.0, 0.0),
                       Complex(1.0, 0.0)};
    AssignFromTo(mat1, mat3, 40);
    calc::MatrixMultiply(mat3, mat4, 10, 1);

    sput_fail_unless(Equal(mat3[0], Complex(7.0, 0.0)), "Check: matrix multiply");
    sput_fail_unless(Equal(mat3[5], Complex(7.0, 0.0)), "Check: matrix multiply");
    sput_fail_unless(Equal(mat3[10], Complex(4.0, 0.0)), "Check: matrix multiply");
    sput_fail_unless(Equal(mat3[22], Complex(12.0, 0.0)), "Check: matrix multiply");
}

void TestInverse()
{
    Complex mat1[40];
    for (int i = 0; i < 10; i++) {
        mat1[0 * 10 + i] = Complex(1.0, 0.0);
        mat1[1 * 10 + i] = Complex(2.0, 0.0);
        mat1[2 * 10 + i] = Complex(2.0, 0.0);
        mat1[3 * 10 + i] = Complex(3.0, 0.0);
    }
    Complex mat2[40];
    AssignFromTo(mat1, mat2, 40);
    calc::MatrixInverse(mat2, 10);
    calc::MatrixMultiply(mat2, mat1, 10);

    sput_fail_unless(Equal(mat2[0], Complex(1.0, 0.0)), "Check: matrix multiply");
    sput_fail_unless(Equal(mat2[5], Complex(1.0, 0.0)), "Check: matrix multiply");
    sput_fail_unless(Equal(mat2[10], Complex(0.0, 0.0)), "Check: matrix multiply");
    sput_fail_unless(Equal(mat2[20], Complex(0.0, 0.0)), "Check: matrix multiply");
    sput_fail_unless(Equal(mat2[30], Complex(1.0, 0.0)), "Check: matrix multiply");
}

void TestG()
{
}

void TestW()
{
}
