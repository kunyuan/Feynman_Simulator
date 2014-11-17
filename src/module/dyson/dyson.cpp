//
//  dyson.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/4/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "dyson.h"
#include "module/parameter/parameter.h"
#include "module/weight/weight.h"

using namespace dyson;
using namespace weight;

Complex* Inverse(Complex*);

bool Dyson::BuildNew(para::ParaDyson &para, Weight &weight)
{
    Beta = para.Beta;
    G = weight.G;
    W = weight.W;
    Sigma = weight.Sigma;
    Polar = weight.Polar;
    return true;
}

void Dyson::DeriveG()
{
    //    G->_Weight();
}

/**
 *  2*2 matrix inverse ( 0  1
 *                       2  3)
 *
 *  @param Complex* original matrix
 *
 *  @return  the inversed matrix
 */
Complex* Inverse(Complex* matrix)
{
    Complex denominator = matrix[0]*matrix[3]-matrix[1]*matrix[2];
    Complex tmp;
    tmp = matrix[3]/denominator;
    matrix[3] = matrix[0]/denominator;
    matrix[0] = tmp;
    matrix[1] *= -1.0/denominator;
    matrix[2] *= -1.0/denominator;
    return matrix;
}

