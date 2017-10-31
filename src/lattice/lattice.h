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
 *  class Lattice includes three set of vectors: 1) LatticeVec (unit cell lattice vector); 2)ReLatticeVec (reciprocal lattice vector for k); 3) SublatticeVec (vectors between different sublattices in the same unit cell).
 */
class Lattice {
public:
    int Dimension;
    int Vol;
    int SublatVol;
    Vec<int> Size;

    Lattice(const Vec<int>& size = Vec<int>(4), int NSublat = 2);
    void Initialize(const Vec<int>& size, int NSublat);

    int Vec2Index(const Vec<int>&) const;
    int Vec2Index(std::initializer_list<int> list) const;
    Vec<int> Index2Vec(int) const;
    int CoordiIndex(const Site& in, const Site& out) const;
    void Shift(Vec<int>& vec) const;
};

int TestLattice();

#endif /* defined(__Fermion_Simulator__lattice__) */
