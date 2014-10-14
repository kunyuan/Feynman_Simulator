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

class GLine {
    friend std::ostream &operator<<(std::ostream &os, GLine &r);
    friend std::istream &operator>>(std::istream &is, GLine &r);

  private:
    int SublatticeIndex;
    int SpinIndex;

  public:
    int Name;
    int Vertex[2];
    spin Spin[2];
    Complex Weight;
    std::string PrettyString();
};

class WLine {
    friend std::ostream &operator<<(std::ostream &os, WLine &r);
    friend std::istream &operator>>(std::istream &is, WLine &r);

  private:
    int SublatticeIndex;

  public:
    int Name;
    int Vertex[2];
    spin Spin[2][2];
    Complex Weight;
    std::string PrettyString();
};

class Vertex {
    friend std::ostream &operator<<(std::ostream &os, Vertex &r);
    friend std::istream &operator>>(std::istream &is, Vertex &r);

  public:
    int Name;
    int G[2];
    int W;
    int Sublattice;
    double tau;
    int r; //TODO: use real space coordinate here
    std::string PrettyString();
};

#endif /* defined(__Fermion_Simulator__component__) */
