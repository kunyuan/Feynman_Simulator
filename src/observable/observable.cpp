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

inline int TauToBin_Sym(real tau)
{
    if(tau>0)
}

Complex SigmaWeight::Weight(Distance dR, real dtau, spin SpinIn, spin SpinOut)
{
    Complex _Weight[SPIN2][NSublattice2][Vol][MAX_BIN];
    return _Weight[SpinIndex(SpinIn, SpinOut)][dR.dSublattice][dR.dCoordinate][;
}



