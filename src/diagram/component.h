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
    friend bool operator==(const GLine &, const GLine &);
    friend bool operator!=(const GLine &, const GLine &);

  public:
    name Name;
    vertex nVer[2];
    int K;
    Complex Weight;
};

class WLine {
    friend std::ostream &operator<<(std::ostream &os, WLine &r);
    friend bool operator==(const WLine &, const WLine &);
    friend bool operator!=(const WLine &, const WLine &);

  public:
    name Name;
    bool IsWorm;
    vertex nVer[2];
    int K;
    Complex Weight;
};

class Vertex {
    friend std::ostream &operator<<(std::ostream &os, Vertex &r);
    friend bool operator==(const Vertex &, const Vertex &);
    friend bool operator!=(const Vertex &, const Vertex &);

  public:
    name Name;
    gLine nG[2];
    wLine nW;
    spin Spin[2]; // IN/OUT spins
    Site R;
    double Tau;
};

class WormClass {
  public:
    bool Exist;
    vertex Ira, Masha; //extra line: Ira---"k,dSpin"--->Masha
    int K;
    int dSpin;

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
