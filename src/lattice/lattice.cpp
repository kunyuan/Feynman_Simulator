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
    SublatVol2 = NSublat * NSublat;
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
 *  get a vector within system size L
 *
 *  @param vec the initial vector
 *
 *  @return new variable within [0, L]
 */
Vec<int> Lattice::Shift(const Vec<int>& vec) const
{
    Vec<int> newvec(vec);
    for (int i = 0; i < D; i++) {
        if (vec[i] < 0)
            newvec[i] = vec[i] + Size[i];
        if (vec[i] >= Size[i])
            newvec[i] = vec[i] - Size[i];
    }
    return newvec;
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
 *  get the sublattice index from In to Out
 *
 *  @para Out: the outgoing sublattice In: the incoming sublattice
 *
 *  @return [0, NSublattice**2-1]
 */
int Lattice::Sublat2Index(int InSub, int OutSub) const
{
    return InSub * SublatVol + OutSub;
}

int Lattice::Index2Sublat(int index, int dir) const
{
    if (dir == IN)
        return index / SublatVol;
    else
        return index % SublatVol;
}

std::tuple<int, int> Lattice::Index2Sublat(int index) const
{
    return std::make_tuple(index / SublatVol, index % SublatVol);
}

bool Lattice::IsOnSameSubLat(int Index)
{
    return Index2Sublat(Index, IN) == Index2Sublat(Index, OUT);
}

/**
 *  get the corresponding site of a name [0, Vol)
 *
 *  @param i name of a site [0, NSublattice*Vol)
 *
 *  @return a site struct
 */

Distance Lattice::Dist(const Site& SiteIn, const Site& SiteOut) const
{
    int sub = Sublat2Index(SiteIn.Sublattice, SiteOut.Sublattice);
    int coord = Vec2Index(Shift(SiteOut.Coordinate - SiteIn.Coordinate));
    return Distance(sub, coord);
}

std::tuple<Site, Site> Lattice::GetSite(const Distance& dis) const
{
    return std::make_tuple(Site(Index2Sublat(dis.SublatIndex, IN), Vec<int>(0)),
                           Site(Index2Sublat(dis.SublatIndex, OUT),
                                Index2Vec(dis.CoordiIndex)));
}

Site Lattice::GetSite(const Distance& dis, int direction) const
{
    if (direction == IN)
        return Site(Index2Sublat(dis.SublatIndex, direction), Vec<int>(0));
    else
        return Site(Index2Sublat(dis.SublatIndex, direction), Index2Vec(dis.CoordiIndex));
}
