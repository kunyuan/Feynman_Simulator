//
//  lattice.h
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/6/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Fermion_Simulator__lattice__
#define __Fermion_Simulator__lattice__

#include "vector.h"

/**
 *  class Site defines all the vertexes on the lattice using a vector of the unit cell: Coordinate and the sublattice number: SubLattice
 */
class Site{
public:
    int SubLattice;
    Vec<int> Coordinate;
    
    Site():SubLattice(0){}
    Site(const int& i):SubLattice(i){}
    int GetName();
    Vec<real> GetVec();
    
    bool operator==(const Site& v2)
    {
        if(SubLattice!=v2.SubLattice) return false;
        if(!(Coordinate==v2.Coordinate)) return false;
        return true;
    }
};

Site GetSite(const int&);

/**
 *  class Distance defines Site1-Site2 using the vector: Dr and sublattice of Site1 and Site2: SubLattice[0] and SubLattice[1].
 */

class Distance{
public:
    int SubLattice[2];
    Vec<int> Dr;
    
    Distance(){}
    Distance(const Site& i, const Site& j)
    {
        Dr=j.Coordinate-i.Coordinate;
        SubLattice[0]=i.SubLattice;
        SubLattice[1]=j.SubLattice;
    }
    
    Vec<real> GetVec();
    Distance Mirror();
};

Distance operator-(const Site& S1, const Site& S2);


/**
 *  class Lattice includes three set of vectors: 1) LatticeVec (unit cell lattice vector); 2)ReLatticeVec (reciprocal lattice vector for k); 3) SubLatticeVec (vectors between different sublattices in the same unit cell).
 */
class Lattice{
public:
    Vec<real> LatticeVec[D];
    Vec<real> ReLatticeVec[D];
    Vec<real> SubLatticeVec[NSublattice];
    
    Lattice() {Initialize();}
    void PlotLattice();
    
private:
    void Initialize();
};

extern Lattice lattice;

int Mirror(int, const int&);
int TestLattice();

#endif /* defined(__Fermion_Simulator__lattice__) */
