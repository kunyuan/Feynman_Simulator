//
//  lattice.h
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/6/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Fermion_Simulator__lattice__
#define __Fermion_Simulator__lattice__

#include "utility/vector.h"
#include <tuple>

int GetSublatIndex(int, int);

/**
 *  class Site defines all the vertexes on the lattice using a vector of the unit cell: Coordinate and the sublattice number: Sublattice
 */
class Site {
public:
    int Sublattice;
    Vec<int> Coordinate;

    Site(int sub = 0, Vec<int> vec = Vec<int>())
        : Sublattice(sub)
        , Coordinate(vec)
    {
    }
};

bool operator==(const Site& v1, const Site& v2);
bool operator!=(const Site& v1, const Site& v2);

/**
 *  class Distance defines Site1-Site2 using the vector: Dr and sublattice of Site1 and Site2: Sublattice[0] and Sublattice[1].
 */

class Distance {
public:
    int SublatIndex;
    int CoordiIndex;

    Distance()
        : SublatIndex(0)
        , CoordiIndex(0)
    {
    }

    Distance(int sublat, int coordi)
        : SublatIndex(sublat)
        , CoordiIndex(coordi)
    {
    }

    bool operator==(const Distance& v2)
    {
        if (SublatIndex != v2.SublatIndex)
            return false;
        if (CoordiIndex != v2.CoordiIndex)
            return false;
        return true;
    }

    bool operator!=(const Distance& v2)
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
    int SublatVol2;
    Vec<int> Size;

    Lattice(const Vec<int>& size = Vec<int>(4), int NSublat = 2);
    void Initialize(const Vec<int>& size, int NSublat);

    int Vec2Index(const Vec<int>&) const;
    int Vec2Index(std::initializer_list<int> list) const;
    Vec<int> Index2Vec(int) const;

    int Sublat2Index(int, int) const;
    int Index2Sublat(int, int direction) const;
    /**
    *  return a tuple of Sublattice
    *
    *  @return tuple<Sublattice IN, Sublattice OUT>
    */
    std::tuple<int, int> Index2Sublat(int index) const;
    //return true if the Index represent a pair between the same sublattice
    bool IsOnSameSubLat(int Index);

    /**
    *  return a tuple of Site
    *
    *  @return tuple<Site IN, Site OUT>
    */
    std::tuple<Site, Site> GetSite(const Distance& dis) const;
    Site GetSite(const Distance& dis, int direction) const;

    Distance Dist(const Site&, const Site&) const;

private:
    Vec<int> Shift(const Vec<int>& vec) const;
};

int TestLattice();

#endif /* defined(__Fermion_Simulator__lattice__) */
