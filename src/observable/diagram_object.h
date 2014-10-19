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
#include "array.h"

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
    Array::array4<Complex> * _Weight;
    Array::array4<Complex> * _WeightAccu;
    Array::array4<Complex> * _WeightSquareAccu;
    
public:
    SigmaWeight(int Vol, int Sublattice);
    ~SigmaWeight();
    void UpdateWeight();
    inline Complex Weight(const Distance& dR, real dtau,spin,spin);
    Estimate<Complex> WeightWithError(const Distance& dR, real dtau, spin, spin);
    inline void Measure(const Complex &weight,const Distance& dR, real dtau, spin, spin);
//    std::string PrettyString();
    void SaveState(std::string);
    bool ReadState(std::string);
};

int TestObservable();
#endif /* defined(__Feynman_Simulator__observable__) */
