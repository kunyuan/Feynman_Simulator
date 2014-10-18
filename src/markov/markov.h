//
//  markov.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/17/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__markov__
#define __Feynman_Simulator__markov__

#include "environment.h"

class Markov
{
private:
    EnvMoneCarlo *_Env;
    Diagram *_Diag;
    RandomFactory *_RNG;
    
public:
    Markov(EnvMoneCarlo *);
    void Hop(const int&);
    void CreateWorm();
    void DeleteWorm();
};

#endif /* defined(__Feynman_Simulator__markov__) */
