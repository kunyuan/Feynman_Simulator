//
//  crc32_test.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/23/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "crc32.h"
#include "sput.h"

void Test_CRC();

int TestCRC32()
{
    sput_start_testing();
    sput_enter_suite("Test crc32...");

    sput_run_test(Test_CRC);

    sput_finish_testing();
    return sput_get_return_value();
}

void Test_CRC()
{
    const char *buffer = "123456789";
    uint32_t crc = 0;
    crc = crc32((uint32_t)0L, (unsigned char *)buffer, (size_t)9);
    sput_fail_unless(crc == 0xCBF43926, "crc32 check");

    crc = crc32(0L, (unsigned char *)NULL, (size_t)1);
    sput_fail_unless(crc == 0, "crc32 check");
}