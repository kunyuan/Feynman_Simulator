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
#include <iostream>

const int MAX_BIN=128;
/**
*  Estimate gives a estimation of certain quantity with it's error bar'
*/

//TODO: Add fitting function here

class WeightG
{
private:
    real _Norm;
    Complex _Matrix[MAX_BIN];
    Complex _ErrorMatrix[MAX_BIN];
public:
    WeightG();
    static Complex& WeightOfBareG(int dr, real dt, spin spin_in, spin spin_out);
    Complex& Weight(int dr, real dt, spin spin_in, spin spin_out);
    Estimate<Complex>& WeightWithError(int dr, real dt, spin spin_in, spin spin_out);
    Complex& operator()(int dr, real dt, spin spin_in, spin spin_out);
    void AddStatistics(Complex &weight, real norm);
    void Write2Binary(std::ostream &);
    void Write2Readable(std::ostream &);
    bool ReadFromReadable(std::istream &);
    bool ReadFromBinary(std::istream &);
};

#endif /* defined(__Feynman_Simulator__observable__) */
