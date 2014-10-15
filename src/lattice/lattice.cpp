//
//  lattice.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/6/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//
#include "definition_global.h"
#include "lattice.h"
#include "math.h"
#include "convention.h"
#include "cnpy.h"
#include <iostream>

Lattice lattice;

/**
 *  initialize the lattice vectors
 */
void Lattice::Initialize()
{
    //Square Lattice
    LatticeVec[0][0]=1.0;
    LatticeVec[0][1]=0.0;
    
    LatticeVec[1][0]=0.0;
    LatticeVec[1][1]=1.0;
    
    ReLatticeVec[0][0]=2.0*Pi;
    ReLatticeVec[0][1]=0.0;
    ReLatticeVec[1][0]=0.0;
    ReLatticeVec[1][1]=2.0*Pi;
    
    SubLatticeVec[0][0]=0.0;
    SubLatticeVec[0][1]=0.0;
    SubLatticeVec[1][0]=0.5;
    SubLatticeVec[1][1]=0.5;
    
    //Lattice Honeycomb
//    LatticeVec[0][0]=0.0;
//    LatticeVec[0][1]=1.0;
//    
//    LatticeVec[1][0]=sqrt(3.0)/2.0;
//    LatticeVec[1][1]=-0.5;
//    
//    ReLatticeVec[0][0]=2.0*Pi/sqrt(3.0);
//    ReLatticeVec[0][1]=2.0*Pi;
//    ReLatticeVec[1][0]=4.0*Pi/sqrt(3.0);
//    ReLatticeVec[1][1]=0.0;
//    
//    SubLatticeVec[0][0]=0.0;
//    SubLatticeVec[0][1]=0.0;
//    SubLatticeVec[1][0]=1.0/2.0/sqrt(3.0);
//    SubLatticeVec[1][1]=0.5;
    
}

/**
 *  get the name for a site on the lattice
 *
 *  @return Name for the Site in [0, Vol)
 */
int Site::GetName()
{
    int cell, layer;
    cell = 0;
    for(int i=D-1; i>=0; i--)
    {
        layer = Coordinate[i];
        for(int j=i-1; j>=0; j--)
            layer = layer*L[j];
        cell = cell + layer;
    }
    return cell*NSublattice+SubLattice;
}

/**
 *  get the real vector for each site
 *
 *  @return a vector which defines the site's coordinate on the lattice
 */
Vec<real> Site::GetVec()
{
    Vec<real> vec(0.0);
    for(int i=0; i<D; i++)
        vec += lattice.LatticeVec[i]*Coordinate[i];
    vec += lattice.SubLatticeVec[SubLattice];
    return vec;
}

Site GetSite(const int& i)
{
    int name, layer;
    Site s;
    name = i/2;
    s.SubLattice = i-name*2;
    
    for(int j=D-1; j>=0; j--)
    {
        layer = 1;
        for(int k=j-1; k>=0; k--)
            layer = layer*L[k];
        s.Coordinate[j] = name/layer;
        name = name - s.Coordinate[j]*layer;
    }
    return s;
}

/**
 *  get the real vector for the distance between two sites
 *
 *  @return a real vector $\vec{r_2}-\vec{r_1}$
 */
Vec<real> Distance::GetVec()
{
    Vec<real> vec(0.0);
    for(int i=0; i<D; i++)
        vec += lattice.LatticeVec[i]*Dr[i];
    vec += lattice.SubLatticeVec[SubLattice[1]]-lattice.SubLatticeVec[SubLattice[0]];
    return vec;
}

/**
 *  get a new distance which is symmetric with current distance and within the first quarter of the lattice
 *
 *  @return a distance within [0, (L/2)]
 */
Distance Distance::Mirror()
{
    Distance dis;
    dis=*this;
    for(int i=0; i<D; i++)
        dis.Dr[i] = ::Mirror(Dr[i], L[i]);
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
int Mirror(int i, const int& L)
{
    if(i<0) i=abs(i);
    if(i>L/2) i = L-i;
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
Distance operator-(const Site& i, const Site& j)
{
    Distance dis;
    dis.Dr=i.Coordinate-j.Coordinate;
    dis.SubLattice[1] = i.SubLattice;
    dis.SubLattice[0] = j.SubLattice;
    return dis;
}

void Lattice::PlotLattice()
{
    const unsigned int N=Vol;
    //save it to file
    const unsigned int shape[] = {N, (unsigned int)D};
    real data[N*D];
    Site s;
    for(int i; i<N; i++)
    {
        s=GetSite(i);
        for(int j=0; j<D; j++)
            data[i*D+j]=s.GetVec()[j];
        std::cout<< data[i*D]<<' '<<data[i*D+1]<<std::endl;
    }
    cnpy::npy_save("sq.npy",data,shape,2,"w");
}
