//
//  lattice.h
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/6/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Fermion_Simulator__lattice__
#define __Fermion_Simulator__lattice__

#include "vector.h"

int GetSublatIndex(int, int);

/**
 *  class Site defines all the vertexes on the lattice using a vector of the unit cell: Coordinate and the sublattice number: Sublattice
 */
class Site {
  public:
    int Sublattice;
    Vec<int> Coordinate;

    Site()
        : Sublattice(0), Coordinate(0)
    {
    }

    Site(int sub, Vec<int> vec)
        : Sublattice(sub), Coordinate(vec)
    {
    }

    bool operator==(const Site &v2)
    {
        if (Sublattice != v2.Sublattice)
            return false;
        if (Coordinate != v2.Coordinate)
            return false;
        return true;
    }

    bool operator!=(const Site &v2)
    {
        if (*this == v2)
            return false;
        return true;
    }
};

/**
 *  class Distance defines Site1-Site2 using the vector: Dr and sublattice of Site1 and Site2: Sublattice[0] and Sublattice[1].
 */

class Distance {
  public:
    int SublatIndex;
    int CoordiIndex;

    Distance()
    {
    }

    Distance(int sublat, int coordi)
        : SublatIndex(sublat), CoordiIndex(coordi)
    {
    }

    int Sublattice(const int &) const;

    bool operator==(const Distance &v2)
    {
        if (SublatIndex != v2.SublatIndex)
            return false;
        if (CoordiIndex != v2.CoordiIndex)
            return false;
        return true;
    }

    bool operator!=(const Distance &v2)
    {
        if (*this == v2)
            return false;
        return true;
    }
};

/**
 *  class Lattice includes three set of vectors: 1) LatticeVec (unit cell lattice vector); 2)ReLatticeVec (reciprocal lattice vector for k); 3) SublatticeVec (vectors between different sublattices in the same unit cell).
 */
class Lattice {
  public:
    int Dimension;
    int Vol;
    int SublatVol;
    Vec<int> Size;
    Vec<real> LatticeVec[D];
    Vec<real> ReciprocalLatticeVec[D];
    Vec<real> SublatticeVec[NSublattice];

    Lattice();
    void PlotLattice();

    Vec<int> &Shift(Vec<int> &);

    int ToIndex(const Vec<int> &);
    Vec<int> ToVec(int);

    Site GetSite(int name);
    int GetName(const Site &);
    int ToIndex(const Site &);
    Vec<real> GetRealVec(const Site &);

    Distance Distance(const Site &, const Site &);
    Vec<int> Coordinate(const class Distance &);
    Vec<real> GetRealVec(const class Distance &);

  private:
    void Initialize();
};

int TestLattice();

#endif /* defined(__Fermion_Simulator__lattice__) */
