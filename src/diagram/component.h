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
class WormClass;

typedef Vertex *vertex;
typedef GLine *gLine;
typedef WLine *wLine;
typedef int name;

class GLine {
    friend std::ostream &operator<<(std::ostream &os, GLine &r);

  public:
    name Name;
    bool IsMeasure;
    vertex nVer[2];
    int K;
    Complex Weight;
};

class WLine {
    friend std::ostream &operator<<(std::ostream &os, WLine &r);

  public:
    name Name;
    bool IsWorm;
    bool IsDelta;
    bool IsMeasure;
    vertex nVer[2];
    int K;
    Complex Weight;
};

class Vertex {
    friend std::ostream &operator<<(std::ostream &os, Vertex &r);

  public:
    name Name;
    gLine nG[2];
    wLine nW;
    spin Spin[2]; // IN/OUT spins
    Site R;
    double Tau;
    int Dir;
};

class WormClass {
  public:
    bool Exist;
    vertex Ira, Masha; //extra line: Ira---"k,dSpin"--->Masha
    int K;
    int dSpin;
    real Weight;

    WormClass()
        : Exist(false), Ira(nullptr), Masha(nullptr), K(0), dSpin(0)
    {
    }
    WormClass(vertex ira, vertex masha, int dk, int s)
        : Exist(true), Ira(ira), Masha(masha), K(dk), dSpin(s)
    {
    }
};

#endif /* defined(__Fermion_Simulator__component__) */
