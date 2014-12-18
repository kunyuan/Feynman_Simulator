//
//  serialization.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 12/17/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__serialization__
#define __Feynman_Simulator__serialization__

#include <Python/Python.h>
//#include "numpy/arrayobject.h"

#include <string>

class Serialization {
public:
    Serialization();
    ~Serialization();
    void Set(std::string key, int value);
    int Get(std::string key);
    void Print();

private:
    PyObject* _Dict;
};

int TestSerialization();

#endif /* defined(__Feynman_Simulator__serialization__) */
