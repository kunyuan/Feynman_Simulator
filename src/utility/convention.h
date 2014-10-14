//
//  convention.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/13/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef Feynman_Simulator_convention_h
#define Feynman_Simulator_convention_h

#ifndef real
#define real double
#endif

enum spin { DOWN,
            UP };
//Spin DOWN: 0,  Spin UP:1
const int IN = 0;
const int OUT = 1;

//lattice information
const int D=2;
const int NSublattice=2;
const int L[D]={4,4};
const int Vol=L[0]*L[1];

#endif
