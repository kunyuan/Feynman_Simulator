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
#define FlipSpin(x) (1-x)
//Spin DOWN: 0,  Spin UP:1
const int SPIN=2;
const int SPIN2=SPIN*SPIN;
const int SPIN3=SPIN2*SPIN;
const int SPIN4=SPIN2*SPIN2;

const int IN = 0;
const int OUT = 1;
#define FlipDir(x) (1-x)

//lattice information
const int D=2;
const int NSublattice=2;
const int NSublattice2=NSublattice*NSublattice;
const int Lx=16;
const int Ly=32;
const int L[D]={Lx,Ly};
const int Vol=Lx*Ly;

const real Pi=3.14159265359;

#endif
