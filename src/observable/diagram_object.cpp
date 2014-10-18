//
//  observable.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/13/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "diagram_object.h"

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
    return 0;
}

Complex SigmaWeight::Weight(Distance dR, real dtau, spin SpinIn, spin SpinOut)
{
    return _Weight[SpinIndex(SpinIn, SpinOut)][dR.dSublattice][0][TauToBin_Sym(dtau)];
    //TODO: need dR.dCoordinate as a lattice Index here
}



