//
//  convention.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/13/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef Feynman_Simulator_convention_h
#define Feynman_Simulator_convention_h

typedef double real;
//#define NDEBUG
//turn off all assert

const bool DEBUGMODE = true;
//#define NDEBUG
//define NDEBUG will turn off debug checking, including the boundary check in array.h

const real PI = 3.1415926535897932384626433832795;

enum spin { DOWN,
            UP };

#define FLIP(x) spin(1 - x)

//Spin DOWN: 0,  Spin UP:1
const int SPIN = 2;
const int SPIN2 = SPIN * SPIN;
const int SPIN3 = SPIN2 * SPIN;
const int SPIN4 = SPIN2 * SPIN2;

const int IN = 0;
const int OUT = 1;
#define INVERSE(x) (1 - x)

const int MAX_ORDER = 15;

//define the model you want to simualte here
enum model { TEST,
             DIAGRAMCOUNTER,
             J1J2,
             HUBBARD };
const model MODEL = J1J2;

//define your lattice here

enum lattice {
    SQUARE,
    CHECKBOARD,
    HONEYCOMB,
    SIMPLE_CUBIC
};
#define LATTICE_INDEX 1

#if LATTICE_INDEX == 0
const lattice LATTICE = SQUARE;
const int D = 2;
const int NSublattice = 1;
#elif LATTICE_INDEX == 1
const lattice LATTICE = CHECKBOARD;
const int D = 2;
const int NSublattice = 2;
#elif LATTICE_INDEX == 2
const lattice LATTICE = HONEYCOMB;
const int D = 2;
const int NSublattice = 2;
#elif LATTICE_INDEX == 3
const lattice LATTICE = SIMPLE_CUBIC;
const int D = 3;
const int NSublattice = 1;
#endif

const int NSublattice2 = NSublattice * NSublattice;

#endif
