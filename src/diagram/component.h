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
  public:
    name Name;
    bool IsMeasure;
    vertex nVer[2];
    int K;
    Complex Weight;
    spin Spin();
    spin Spin(int dir);
    void FlipSpin();
    int Sublattice(int dir);
    vertex NeighVer(int dir);
    std::string PrettyString();
};

class WLine {
  public:
    name Name;
    bool IsWorm;
    bool IsDelta;
    bool IsMeasure;
    vertex nVer[2];
    int K;
    Complex Weight;
    spin Spin(int, int);
    void FlipSpin();
    int Sublattice(int dir);
    vertex NeighVer(int dir);
    std::string PrettyString();
};

class Vertex {
  public:
    name Name;
    gLine nG[2];
    wLine nW;
    Site R;
    spin _spin[2]; // IN/OUT spins
    double Tau;
    int Dir;
    spin Spin(int);
    void SetSpin(spin *);
    spin *Spin();
    int Sublattice();
    gLine NeighG(int dir);
    wLine NeighW();
    std::string PrettyString();
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
