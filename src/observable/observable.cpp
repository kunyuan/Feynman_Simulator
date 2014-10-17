//
//  observable.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/13/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "observable.h"

inline int SpinIndex(spin SpinIn, spin SpinOut)
{
    return SpinIn*SPIN+SpinOut;
}

inline int SpinIndex(spin* TwoSpinIn, spin* TwoSpinOut)
{
    return TwoSpinIn[0]*SPIN3+TwoSpinIn[1]*SPIN2+
        TwoSpinOut[0]*SPIN+TwoSpinOut[1];
}

Complex SigmaWeight::Weight(Distance dR, real dtau, spin SpinIn, spin SpinOut)
{
    return Complex(0.0,0.0);
}



