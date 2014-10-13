//
//  observable.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/13/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__observable__
#define __Feynman_Simulator__observable__

#include "complex.h"

class WeightG
{
    Complex& Weight();
    Complex& operator()();
};

#endif /* defined(__Feynman_Simulator__observable__) */
