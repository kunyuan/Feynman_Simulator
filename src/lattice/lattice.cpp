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
//    LatticeVec[1][1]=0.5;
//    
//    ReLatticeVec[0][0]=2.0*Pi/sqrt(3.0);
//    ReLatticeVec[0][1]=-2.0*Pi;
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
