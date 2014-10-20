//
//  lattice.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/6/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//
#include "lattice.h"
#include "math.h"
#include "convention.h"
#include "cnpy.h"
#include "abort.h"
#include <iostream>

Lattice::Lattice()
{
    Dimension = D;
    Vol=1;
    for (int i = 0; i < D; i++)
    {
        Size[i] = L[i];
        Vol*=L[i];
    }
    SublatVol=NSublattice;
    Initialize();
}
/**
 *  initialize the lattice vectors
 */
void Lattice::Initialize()
{
    //Square Lattice
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
 *  get the name for a site on the lattice
 *
 *  @return Name for the Site in [0, NSublattice*Vol)
 */
int Site::GetName()
{
    return Coordinate.ToIndex() * NSublattice + Sublattice;
}

/**
 *  get the corresponding site of a name [0, Vol)
 *
 *  @param i name of a site [0, NSublattice*Vol)
 *
 *  @return a site struct
 */
Site::Site(const int name)
    : Coordinate(name / 2), Sublattice(name % 2)
{
}

/**
 *  get the sublattice index from In to Out
 *
 *  @para Out: the outgoing sublattice In: the incoming sublattice
 *
 *  @return [0, NSublattice**2-1]
 */
int GetSublatIndex(int In, int Out)
{
    if (DEBUGMODE) {
        if (In < 0 || In >= NSublattice)
            ABORT("Wrong sublattice number!");
        if (Out < 0 || Out >= NSublattice)
            ABORT("Wrong sublattice number!");
    }
    return Out * NSublattice + In;
}

int Distance::Sublattice(const int &dir) const
{
    if (DEBUGMODE && (_SublatIndex < 0 || _SublatIndex >= NSublattice2))
        ABORT("Wrong sublattice index number!");
    if (dir == IN)
        return _SublatIndex % NSublattice;
    else
        return _SublatIndex / NSublattice;
}

Vec<int> Distance::Coordinate() const
{
    if (DEBUGMODE && (_CoordiIndex < 0 || _CoordiIndex >= Vol))
        ABORT("Wrong Coordinate index number!");
    return Vec<int>(_CoordiIndex);
}

/**
 *  get a new distance which is symmetric with current distance and within the first quarter of the lattice
 *
 *  @return a distance within [0, (L/2)]
 */
Distance Distance::Mirror()
{
    Distance dis;
    dis = *this;
    Vec<int> v;
    for (int i = 0; i < D; i++)
        v[i] = ::Mirror(Coordinate()[i], L[i]);
    dis._CoordiIndex = v.ToIndex();
    return dis;
}
/**
 *  get a int with is a symmetric with i within system size L
 *
 *  @param i the initial int
 *  @param L system size
 *
 *  @return new variable within [0, L/2] which is symmetric with i
 */
int Mirror(int i, const int &L)
{
    if (i < 0)
        i = abs(i);
    if (i > L / 2)
        i = L - i;
    return i;
}

/**
 *  Distance = Site1 - Site2
 *
 *  @param i Site1
 *  @param j Site2
 *
 *  @return a new Distance(Site1-Site2)
 */
Distance operator-(const Site &i, const Site &j)
{
    Distance dis;
    dis._SublatIndex = GetSublatIndex(j.Sublattice, i.Sublattice);
    dis._CoordiIndex = (Vec<int>(i.Coordinate - j.Coordinate)).ToIndex();
    return dis;
}

/**
 *  get the real vector for each site
 *
 *  @return a vector which defines the site's coordinate on the lattice
 */
Vec<real> Lattice::GetRealVec(const Site &site)
{
    Vec<real> vec(0.0);
    for (int i = 0; i < D; i++)
        vec += LatticeVec[i] * site.Coordinate[i];
    vec += SublatticeVec[site.Sublattice];
    return vec;
}

/**
 *  get the real vector for the distance between two sites
 *
 *  @return a real vector $\vec{r_2}-\vec{r_1}$
 */
Vec<real> Lattice::GetRealVec(const Distance &dis)
{
    Vec<real> vec(0.0);
    for (int i = 0; i < D; i++)
        vec += LatticeVec[i] * dis.Coordinate()[i];
    vec += SublatticeVec[dis.Sublattice(OUT)] - SublatticeVec[dis.Sublattice(IN)];
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
        Site s(i);
        for (int j = 0; j < D; j++)
            data[i * D + j] = GetRealVec(s)[j];
    }
    cnpy::npy_save("sq.npy", data, shape, 2, "w");
}
