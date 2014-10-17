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
#include "convention.h"
#include "estimate.h"
#include "lattice.h"
#include <iostream>

const int MAX_BIN=128;
/**
*  Estimate gives a estimation of certain quantity with it's error bar'
*/

//TODO: Add fitting function here
//TODO: Add reweighting
class SigmaWeight
{
private:
    real _Norm;
    Complex _Weight[SPIN2][NSublattice2][Vol][MAX_BIN];
    Complex _WeightAccu[SPIN2][NSublattice2][Vol][MAX_BIN];
    Complex _ErrorSquare[SPIN2][NSublattice2][Vol][MAX_BIN];
public:
    void Update();
    Complex Weight(Distance dR, real dtau,spin,spin);
    Estimate<Complex> WeightWithError(Distance dR, real dtau, spin, spin);
    void AddStatistics(Complex &weight, real norm);
//    std::string PrettyString();
    void SaveState(std::string);
    bool ReadState(std::string);
};

int TestObservable();
#endif /* defined(__Feynman_Simulator__observable__) */
