//
//  Estimate.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/14/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__Estimate__
#define __Feynman_Simulator__Estimate__

#include <iostream>
#include "component.h"
/**
*  \brief estimate with mean value and standard error
*/
template <typename T>
class Estimate
{
public:
    Estimate();
    Estimate(T &m, T &e);
    T Mean;
    T Error;
    friend std::ostream& operator<<(std::ostream &, Estimate &);
};

/**
*  \brief maintain a history of observables, where an Estimate object can be calculated
*/

class EstimateKeeper
{
    
};

#endif /* defined(__Feynman_Simulator__Estimate__) */
