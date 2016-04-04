//
//  lattice.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/6/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//
#include "lattice.h"
#include <math.h>
#include <iostream>
#include "../utility/convention.h"
#include "../utility/abort.h"
#include "../utility/utility.h"
Lattice::Lattice(const Vec<int>& size, int NSublat)
{
    Initialize(size, NSublat);
}

void Lattice::Initialize(const Vec<int>& size, int NSublat)
{
    Dimension = D;
    Vol = 1;
    Size = size;
    for (int i = 0; i < D; i++) {
        Vol *= Size[i];
    }
    SublatVol = NSublat;
}

bool operator==(const Site& v1, const Site& v2)
{
    if (v1.Sublattice != v2.Sublattice)
        return false;
    if (v1.Coordinate != v2.Coordinate)
        return false;
    return true;
}

bool operator!=(const Site& v1, const Site& v2)
{
    return !(v1 == v2);
}

/**
 * return vec[0]*L1*L2+vec[1]*L2+vec[2]
 * */
int Lattice::Vec2Index(const Vec<int>& vec) const
{
    int Index = vec[0];
    for (int i = 1; i < D; i++) {
        Index = Index * Size[i] + vec[i];
    }
    return Index;
}

int Lattice::Vec2Index(std::initializer_list<int> list) const
{
    return Vec2Index(Vec<int>(list));
}

Vec<int> Lattice::Index2Vec(int index) const
{
    Vec<int> v(0);
    for (int i = D - 1; i > 0; i--) {
        v[i] = index % Size[i];
        index = index / Size[i];
    }
    v[0] = index;
    return v;
}
/**
 *  get a vector within system size L
 *
 *  @param vec the initial vector
 *
 *  @return new variable within [0, L]
 */
void Lattice::Shift(Vec<int>& vec) const
{
    for (int i = 0; i < D; i++) {
        if (vec[i] < 0)
            vec[i] += Size[i];
        if (vec[i] >= Size[i])
            vec[i] -= Size[i];
    }
}
int Lattice::CoordiIndex(const Site& in, const Site& out) const
{
    auto v = out.Coordinate - in.Coordinate;
    Shift(v);
    return Vec2Index(v);
}
