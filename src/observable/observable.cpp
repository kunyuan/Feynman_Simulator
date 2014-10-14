//
//  observable.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/13/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "observable.h"

Estimator::Estimator():Mean(),Error(){}

/**
*  Estimator constructor
*
*  @param mean __reference__ of a Complex number as the mean value
*  @param error __reference__ of a Complex number as the error
*
*/
Estimator::Estimator(Complex& mean, Complex& error):Mean(mean),Error(error){}


