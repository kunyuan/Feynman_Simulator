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
    if (DEBUGMODE && (SublatIndex < 0 || SublatIndex >= NSublattice2))
        ABORT("Wrong sublattice index number!");
    if (dir == IN)
        return SublatIndex % NSublattice;
    else
        return SublatIndex / NSublattice;
}

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
Vec<int>& Lattice::Shift(Vec<int>& vec )
{
    for(int i=0; i<Dimension; i++)
    {
        if (vec[i] < 0)
            vec[i] = vec[i]+Size[i];
    }
    return vec;
}


int Lattice::ToIndex(const Vec<int>& vec)
{
    int Index = vec[Dimension - 1];
    for (int i = Dimension - 2; i >= 0; i--) {
        Index = Index * Size[i] + vec[i];
    }
    return Index;
}

Vec<int> Lattice::ToVec(int index)
{
    Vec<int> v(0);
    for (int i = 0; i < Dimension - 1; i++) {
        v[i] = index % Size[i];
        index = index / Size[i];
    }
    v[Dimension - 1] = index;
    return v;
}
/**
 *  get the corresponding site of a name [0, Vol)
 *
 *  @param i name of a site [0, NSublattice*Vol)
 *
 *  @return a site struct
 */
Site Lattice::GetSite(int name)
{
    return Site(name%2, ToVec(name/2));
}

/**
 *  get the name for a site on the lattice
 *
 *  @return Name for the Site in [0, NSublattice*Vol)
 */
int Lattice::GetName(const Site& site)
{
    return ToIndex(site) * NSublattice + site.Sublattice;
}

int Lattice::ToIndex(const Site& site)
{
    return ToIndex(site.Coordinate);
}

/**
 *  get the real vector for each site
 *
 *  @return a vector which defines the site's coordinate on the lattice
 */
Vec<real> Lattice::GetRealVec(const Site &site)
{
    Vec<real> vec(0.0);
    for (int i = 0; i < Dimension; i++)
        vec += LatticeVec[i] * site.Coordinate[i];
    vec += SublatticeVec[site.Sublattice];
    return vec;
}

Distance Lattice::Distance(const Site& SiteIn, const Site& SiteOut)
{
    Vec<int> vec = SiteOut.Coordinate-SiteIn.Coordinate;
    return ::Distance(GetSublatIndex((SiteIn.Sublattice), (SiteOut.Sublattice)),
                      ToIndex(Shift(vec)));
}


Vec<int> Lattice::Coordinate(const class Distance& dis)
{
    if (DEBUGMODE && (dis.CoordiIndex < 0 || dis.CoordiIndex >= Vol))
        ABORT("Wrong Coordinate index number!");
    return ToVec(dis.CoordiIndex);
}

/**
 *  get the real vector for the distance between two sites
 *
 *  @return a real vector $\vec{r_2}-\vec{r_1}$
 */
Vec<real> Lattice::GetRealVec(const class Distance &dis)
{
    Vec<real> vec(0.0);
    for (int i = 0; i < Dimension; i++)
        vec += LatticeVec[i] * Coordinate(dis)[i];
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
    const unsigned int shape[] = {N, (unsigned int)Dimension};
    real data[N * Dimension];
    for (int i = 0; i < N; i++) {
        Site s=GetSite(i);
        for (int j = 0; j < Dimension; j++)
            data[i * Dimension + j] = GetRealVec(s)[j];
    }
    cnpy::npy_save("sq.npy", data, shape, 2, "w");
}
