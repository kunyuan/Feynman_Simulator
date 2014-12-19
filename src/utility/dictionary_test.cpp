//
//  serialization_test.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 12/17/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//
#include "dictionary.h"
#include <initializer_list>
#include <iostream>

using namespace std;

int TestDictionary()
{
    Dictionary Port;
    Port.Set("hello", 1);
    cout << Port.Get<int>("hello") << endl;
    Port.LoadByEval("{'a':1,'b':2}");
    vector<int> v = { 1, 1, 1 };
    Port.Set("Vec", v);
    Port.Print();
    return 0;
}
