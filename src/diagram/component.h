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

class Vertex;
class WLine;
class GLine;
class WormClass;

typedef Vertex* vertex;
typedef GLine* gLine;
typedef WLine* wLine;

class GLine {
    friend std::ostream &operator<<(std::ostream &os, GLine &r);
    friend bool operator==(const GLine&, const GLine&);
    friend bool operator!=(const GLine&, const GLine&);

  public:
    int Name;
    Vertex * nVer[2];
    int K;
    Complex Weight;
    std::istream &LoadConfig(std::istream &os);
    std::ostream &SaveConfig(std::ostream &is);
};

class WLine {
    friend std::ostream &operator<<(std::ostream &os, WLine &r);
    friend bool operator==(const WLine&, const WLine&);
    friend bool operator!=(const WLine&, const WLine&);

  public:
    int Name;
    bool IsWorm;
    Vertex* nVer[2];
    int K;
    Complex Weight;
    std::istream &LoadConfig(std::istream &os);
    std::ostream &SaveConfig(std::ostream &is);
};

class Vertex {
    friend std::ostream &operator<<(std::ostream &os, Vertex &r);
    friend bool operator==(const Vertex&, const Vertex&);
    friend bool operator!=(const Vertex&, const Vertex&);

  public:
    int Name;
    GLine* nG[2];
    WLine* nW;
    spin Spin[2]; // IN/OUT spins
    Site R;
    double Tau;
    std::istream &LoadConfig(std::istream &os);
    std::ostream &SaveConfig(std::ostream &is);
};

class WormClass {
  public:
    bool Exist;
    Vertex *Ira, *Masha; //extra line: Ira---"k,dSpin"--->Masha
    int K;
    int dSpin;
    std::istream &LoadConfig(std::istream &os);
    std::ostream &SaveConfig(std::ostream &is);

    WormClass()
        : Exist(false), Ira(nullptr), Masha(nullptr), K(0), dSpin(0)
    {
    }
    WormClass(Vertex* ira, Vertex* masha, int dk, int s)
        : Exist(true), Ira(ira), Masha(masha), K(dk), dSpin(s)
    {
    }
};

#endif /* defined(__Fermion_Simulator__component__) */
