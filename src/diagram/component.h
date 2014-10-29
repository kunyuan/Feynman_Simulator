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

class Worm {
  public:
    bool Exist;
    int Ira, Masha; //extra line: Ira---"k,dSpin"--->Masha
    int K;
    int dSpin;
    std::istream &LoadConfig(std::istream &os);
    std::ostream &SaveConfig(std::ostream &is);

    Worm()
        : Exist(false), Ira(0), Masha(0), K(0), dSpin(0)
    {
    }
    Worm(int ira, int masha, int dk, int s)
        : Exist(true), Ira(ira), Masha(masha), K(dk), dSpin(s)
    {
    }
};

class GLine {
    friend std::ostream &operator<<(std::ostream &os, GLine &r);

  public:
    int Name;
    int Vertex[2];
    int K;
    Complex Weight;
    std::istream &LoadConfig(std::istream &os);
    std::ostream &SaveConfig(std::ostream &is);
};

class WLine {
    friend std::ostream &operator<<(std::ostream &os, WLine &r);

  public:
    int Name;
    bool IsWorm;
    int Vertex[2];
    int K;
    Complex Weight;
    std::istream &LoadConfig(std::istream &os);
    std::ostream &SaveConfig(std::ostream &is);
};

class Vertex {
    friend std::ostream &operator<<(std::ostream &os, Vertex &r);

  public:
    int Name;
    int G[2];
    int W;
    spin Spin[2]; // IN/OUT spins
    Site R;
    double Tau;
    std::istream &LoadConfig(std::istream &os);
    std::ostream &SaveConfig(std::ostream &is);
};

#endif /* defined(__Fermion_Simulator__component__) */
