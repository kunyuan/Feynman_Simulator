//
//  Estimate.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/14/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "estimate.h"
using namespace std;

template class Estimate<Complex>;
template class Estimate<real>;

template <typename T>
Estimate<T>::Estimate():Mean(),Error(){}

/**
*  Estimate constructor of type T
*
*  @param mean __reference__ of a T number as the mean value
*  @param error __reference__ of a T number as the error
*
*/
template <typename T>
Estimate<T>::Estimate(T& mean, T& error):Mean(mean),Error(error){}

ostream& operator<<(ostream &os, Estimate<Complex> & e)
{
    return os;
}