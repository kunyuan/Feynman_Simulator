//
//  component.h
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/7/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Fermion_Simulator__component__
#define __Fermion_Simulator__component__

#include <string>
#include "complex.h"
#include "utility.h"
#include "logger.h"
#include "convention.h"
#include "lattice.h"

class GLine {
    friend std::ostream &operator<<(std::ostream &os, GLine &r);
    friend std::istream &operator>>(std::istream &is, GLine &r);

  public:
    int Name;
    int Vertex[2];
    int K;
    Complex Weight;
};

class WLine {
    friend std::ostream &operator<<(std::ostream &os, WLine &r);
    friend std::istream &operator>>(std::istream &is, WLine &r);

  public:
    int Name;
    bool IsWorm;
    int Vertex[2];
    int K;
    Complex Weight;
};

class Vertex {
    friend std::ostream &operator<<(std::ostream &os, Vertex &r);
    friend std::istream &operator>>(std::istream &is, Vertex &r);

  public:
    int Name;
    int G[2];
    int W;
    spin Spin[2];// IN/OUT spins
    Site R;
    double Tau;
};

#endif /* defined(__Fermion_Simulator__component__) */
