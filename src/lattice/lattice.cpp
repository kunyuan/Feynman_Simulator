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
Lattice::Lattice()
{
    Initialize(Vec<int>(4), LATTICE);
}
Lattice::Lattice(const Vec<int> &size, lattice _Lattice)
{
    Initialize(size, _Lattice);
}

void Lattice::Initialize(const Vec<int> &size, lattice _Lattice)
{
    Dimension = D;
    Vol = 1;
    Size = size;
    for (int i = 0; i < D; i++) {
        Vol *= Size[i];
    }
    SublatVol = NSublattice;
    SublatVol2 = NSublattice * NSublattice;
    switch (_Lattice) {
        case SQUARE:
            ABORT("square lattice has not been implemented yet!");
        case CHECKBOARD:
            _Checkboard();
            break;
        case HONEYCOMB:
            _Honeycomb();
            break;
        case SIMPLE_CUBIC:
            ABORT("simple cubic lattice has not been implemented yet!");
    }
}

bool operator==(const Site &v1, const Site &v2)
{
    if (v1.Sublattice != v2.Sublattice)
        return false;
    if (v1.Coordinate != v2.Coordinate)
        return false;
    return true;
}

bool operator!=(const Site &v1, const Site &v2)
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
Vec<int> Lattice::Shift(const Vec<int> &vec) const
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

int Lattice::Vec2Index(const Vec<int> &vec) const
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
int Lattice::Sublat2Index(int In, int Out) const
{
    return Out * NSublattice + In;
}

int Lattice::Index2Sublat(int index, int dir) const
{
    if (dir == IN)
        return index % NSublattice;
    else
        return index / NSublattice;
}

bool Lattice::IsOnSameSubLat(int Index)
{
    return Index2Sublat(Index, IN) ==
           Index2Sublat(Index, OUT);
}

/**
 *  get the corresponding site of a name [0, Vol)
 *
 *  @param i name of a site [0, NSublattice*Vol)
 *
 *  @return a site struct
 */
Site Lattice::GetSite(int name) const
{
    return Site(name % 2, Index2Vec(name / 2));
}

/**
 *  get the name for a site on the lattice
 *
 *  @return Name for the Site in [0, NSublattice*Vol)
 */
int Lattice::GetName(const Site &site) const
{
    return Vec2Index(site.Coordinate) * NSublattice + site.Sublattice;
}

/**
 *  get the real vector for each site
 *
 *  @return a vector which defines the site's coordinate on the lattice
 */
Vec<real> Lattice::GetRealVec(const Site &site, Vec<int> offset) const
{
    Vec<real> vec(0.0);
    Vec<int> coordinate = site.Coordinate;
    coordinate = Shift(coordinate + offset);
    for (int i = 0; i < D; i++)
        vec += LatticeVec[i] * coordinate[i];
    vec += SublatticeVec[site.Sublattice];
    return vec;
}

Distance Lattice::Dist(const Site &SiteIn, const Site &SiteOut) const
{
    return Distance(Sublat2Index(SiteIn.Sublattice, SiteOut.Sublattice),
                    Vec2Index(Shift(SiteOut.Coordinate - SiteIn.Coordinate)));
}

Site Lattice::GetSite(const Distance &dis, int direction) const
{
    if (direction == IN)
        return Site(Index2Sublat(dis.SublatIndex, direction), Vec<int>(0));
    else
        return Site(Index2Sublat(dis.SublatIndex, direction), Index2Vec(dis.CoordiIndex));
}

/**
 *  get the real vector for the distance between two sites
 *
 *  @return a real vector $\vec{r_2}-\vec{r_1}$
 */
Vec<real> Lattice::GetRealVec(const class Distance &dis, Vec<int> offset) const
{
    Vec<real> vec(0.0);
    Site site_out = GetSite(dis, OUT);
    Site site_in = GetSite(dis, IN);
    Vec<int> coordinate = Shift(site_out.Coordinate + offset);
    for (int i = 0; i < D; i++)
        vec += LatticeVec[i] * coordinate[i];
    vec += SublatticeVec[site_out.Sublattice] - SublatticeVec[site_in.Sublattice];
    return vec;
}

/**
 *  save all the coordinates of the sites on lattice to a plotable python file
 */
void Lattice::PlotLattice()
{
    ofstream os("lattice.py", ios::out);
    const unsigned int N = NSublattice * Vol;
    os << "#[Real Vec tuple, Int Vec tuple, Sublattice]" << endl;
    os << "points=[" << endl;
    Vec<int> offset;
    for (int i = 0; i < D; i++)
        offset[i] = Size[i] / 2 - 1;
    for (int i = 0; i < N; i++) {
        Site site = GetSite(i);
        os << "[" << GetRealVec(site, offset).PrettyString() << ","
           << site.Coordinate.PrettyString() << "," << site.Sublattice
           << "]," << endl;
    }
    os << "]";
    os.close();
}

/**
 *  initialize the lattice vectors
 */
void Lattice::_Checkboard()
{
    ASSERT_ALLWAYS(Dimension == 2 && SublatVol == 2,
                   "Checkboard lattice has D=2 and Sublattice=2");
    ASSERT_ALLWAYS(Size[0] > 1 && Size[1] > 1, "System size must be bigger than 1!");
    //Square Lattice with two sublattices
    LatticeVec[0][0] = 1.0;
    LatticeVec[0][1] = 0.0;

    LatticeVec[1][0] = 0.0;
    LatticeVec[1][1] = 1.0;

    ReciprocalLatticeVec[0][0] = 2.0 * PI;
    ReciprocalLatticeVec[0][1] = 0.0;
    ReciprocalLatticeVec[1][0] = 0.0;
    ReciprocalLatticeVec[1][1] = 2.0 * PI;

    SublatticeVec[0][0] = 0.0;
    SublatticeVec[0][1] = 0.0;
    SublatticeVec[1][0] = 0.5;
    SublatticeVec[1][1] = 0.5;
}

void Lattice::_Honeycomb()
{
    //Lattice Honeycomb
    LatticeVec[0][0] = 0.0;
    LatticeVec[0][1] = 1.0;

    LatticeVec[1][0] = sqrt(3.0) / 2.0;
    LatticeVec[1][1] = -0.5;

    ReciprocalLatticeVec[0][0] = 2.0 * PI / sqrt(3.0);
    ReciprocalLatticeVec[0][1] = 2.0 * PI;
    ReciprocalLatticeVec[1][0] = 4.0 * PI / sqrt(3.0);
    ReciprocalLatticeVec[1][1] = 0.0;

    SublatticeVec[0][0] = 0.0;
    SublatticeVec[0][1] = 0.0;
    SublatticeVec[1][0] = 1.0 / 2.0 / sqrt(3.0);
    SublatticeVec[1][1] = 0.5;
}
