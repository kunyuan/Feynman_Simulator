//
//  serialization_test.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 12/17/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//
#include "serialization.h"
#include <iostream>

using namespace std;

int TestSerialization()
{
    Serialization Port;
    Port.Set("hello", 1);
    cout << Port.Get("hello") << endl;
    Port.Print();
    return 0;
}
