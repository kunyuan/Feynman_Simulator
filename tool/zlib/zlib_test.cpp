//
//  zlib_test.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/23/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "zlib.h"
#include "../sput.h"

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
    uLong crc = crc32(0L, (Bytef *)buffer, (uInt)9);
    sput_fail_unless(crc == 0xCBF43926, "crc32 check");

    const char *buffer1 = "0";
    crc = crc32(0L, (Bytef *)buffer1, (uInt)1);
    sput_fail_unless(crc == 0xF4DBDF21, "crc32 check");
}