//
//  component.h
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/7/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Fermion_Simulator__component__
#define __Fermion_Simulator__component__

#include "../utility/complex.h"
#include "../utility/convention.h"
#include "../lattice/lattice.h"

class GLine;
class WLine;
class Vertex;

typedef Vertex *vertex;
typedef WLine *wLine;
typedef GLine *gLine;

class WormClass {
  public:
    bool Exist;
    int Ira, Masha; //extra line: Ira---"k,dSpin"--->Masha
    int K;
    int dSpin;

    WormClass()
        : Exist(false), Ira(0), Masha(0), K(0), dSpin(0)
    {
    }
    WormClass(int ira, int masha, int dk, int s)
        : Exist(true), Ira(ira), Masha(masha), K(dk), dSpin(s)
    {
    }
};

class GLine {
  public:
    int Name;
    Vertex *nVer[2];
    int K;
    Complex Weight;
};

class WLine {
  public:
    int Name;
    bool IsWorm;
    Vertex *nVer[2];
    int K;
    Complex Weight;
};

class Vertex {
  public:
    int Name;
    GLine *G[2];
    WLine *nW;
    spin Spin[2]; // IN/OUT spins
    Site R;
    double Tau;
};

#endif /* defined(__Fermion_Simulator__component__) */
