//
//  array_test.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 1/20/15.
//  Copyright (c) 2015 Kun Chen. All rights reserved.
//

#include "array.h"
#include "utility/utility.h"
#include "sput.h"

void Test_Array_Access();
int TestArray()
{
    sput_start_testing();
    sput_enter_suite("Test Array");
    sput_run_test(Test_Array_Access);
    sput_finish_testing();
    return sput_get_return_value();
}

void Test_Array_Access()
{
    Array<3> array;
    uint shape[3] = { 2, 3, 2 };
    array.Allocate(shape);
    array.Assign({ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 });
    sput_fail_unless(Equal(array.GetShape(), shape, 3), "Check shape");
    sput_fail_unless(Equal(array({ 0, 0, 1 }), 2), "Check element");
    sput_fail_unless(Equal(array({ 1, 1, 1 }), 10), "Check element");
    array({ 1, 1, 0 }) = 4;
    sput_fail_unless(Equal(array({ 1, 1, 0 }), 4), "Check element");
}
