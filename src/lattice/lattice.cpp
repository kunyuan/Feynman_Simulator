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
#include "../utility/cnpy.h"
#include "../utility/abort.h"
Lattice::Lattice(const Vec<int> &size)
{
    Dimension = D;
    Vol = 1;
    Size = size;
    for (int i = 0; i < D; i++) {
        Vol *= Size[i];
    }
    SublatVol = NSublattice;
    Initialize();
}
/**
 *  initialize the lattice vectors
 */
void Lattice::Initialize()
{
    //TODO: Simple Cubic Lattice with two sublattices

    //Square Lattice with two sublattices
    LatticeVec[0][0] = 1.0;
    LatticeVec[0][1] = 0.0;

    LatticeVec[1][0] = 0.0;
    LatticeVec[1][1] = 1.0;

    ReciprocalLatticeVec[0][0] = 2.0 * Pi;
    ReciprocalLatticeVec[0][1] = 0.0;
    ReciprocalLatticeVec[1][0] = 0.0;
    ReciprocalLatticeVec[1][1] = 2.0 * Pi;

    SublatticeVec[0][0] = 0.0;
    SublatticeVec[0][1] = 0.0;
    SublatticeVec[1][0] = 0.5;
    SublatticeVec[1][1] = 0.5;

    //Lattice Honeycomb
    //    LatticeVec[0][0]=0.0;
    //    LatticeVec[0][1]=1.0;
    //
    //    LatticeVec[1][0]=sqrt(3.0)/2.0;
    //    LatticeVec[1][1]=-0.5;
    //
    //    ReciprocalLatticeVec[0][0]=2.0*Pi/sqrt(3.0);
    //    ReciprocalLatticeVec[0][1]=2.0*Pi;
    //    ReciprocalLatticeVec[1][0]=4.0*Pi/sqrt(3.0);
    //    ReciprocalLatticeVec[1][1]=0.0;
    //
    //    SublatticeVec[0][0]=0.0;
    //    SublatticeVec[0][1]=0.0;
    //    SublatticeVec[1][0]=1.0/2.0/sqrt(3.0);
    //    SublatticeVec[1][1]=0.5;
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
    }
    return newvec;
}

int Lattice::Vec2Index(const Vec<int> &vec) const
{
    int Index = vec[D - 1];
    for (int i = D - 2; i >= 0; i--) {
        Index = Index * Size[i] + vec[i];
    }
    return Index;
}

Vec<int> Lattice::Index2Vec(int index) const
{
    Vec<int> v(0);
    for (int i = 0; i < D - 1; i++) {
        v[i] = index % Size[i];
        index = index / Size[i];
    }
    v[D - 1] = index;
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
Vec<real> Lattice::GetRealVec(const Site &site) const
{
    Vec<real> vec(0.0);
    for (int i = 0; i < D; i++)
        vec += LatticeVec[i] * site.Coordinate[i];
    vec += SublatticeVec[site.Sublattice];
    return vec;
}

Distance Lattice::Dist(const Site &SiteIn, const Site &SiteOut) const
{
    return Distance(Sublat2Index(SiteIn.Sublattice, SiteOut.Sublattice),
                    Vec2Index(Shift(SiteOut.Coordinate - SiteIn.Coordinate)));
}

int Lattice::GetSublat(const Distance &dis, int dir) const
{
    return Index2Sublat(dis.SublatIndex, dir);
}

Vec<int> Lattice::GetVec(const class Distance &dis) const
{
    if (DEBUGMODE && (dis.CoordiIndex < 0 || dis.CoordiIndex >= Vol))
        ABORT("Wrong Coordinate index number!");
    return Index2Vec(dis.CoordiIndex);
}

/**
 *  get the real vector for the distance between two sites
 *
 *  @return a real vector $\vec{r_2}-\vec{r_1}$
 */
Vec<real> Lattice::GetRealVec(const class Distance &dis) const
{
    Vec<real> vec(0.0);
    for (int i = 0; i < D; i++)
        vec += LatticeVec[i] * GetVec(dis)[i];
    vec += SublatticeVec[GetSublat(dis, OUT)] - SublatticeVec[GetSublat(dis, IN)];
    return vec;
}

/**
 *  save all the coordinates of the sites on lattice to a plotable python file
 */
void Lattice::PlotLattice()
{
    const unsigned int N = NSublattice * Vol;
    //save it to file
    const unsigned int shape[] = {N, (unsigned int)D};
    real data[N * D];
    for (int i = 0; i < N; i++) {
        Site s = GetSite(i);
        for (int j = 0; j < D; j++)
            data[i * D + j] = GetRealVec(s)[j];
    }
    cnpy::npy_save("sq.npy", data, shape, 2, "w");
}
