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
 *  class Site defines all the vertexes on the lattice using a vector of the unit cell: Coordinate and the sublattice number: Sublattice
 */
class Site{
public:
    int Sublattice;
    Vec<int> Coordinate;
    
    Site():Sublattice(0){}
    Site(const int name);
    int GetName();
    Vec<real> GetVec();
    
    bool operator==(const Site& v2)
    {
        if(Sublattice!=v2.Sublattice) return false;
        if(Coordinate!=v2.Coordinate)) return false;
        return true;
    }
    
    bool operator!=(const Site& v2)
    {
        if(*this==v2) return false;
        return true;
    }
};
int GetSublatticeIndex(const int&, const int&);

/**
 *  class Distance defines Site1-Site2 using the vector: Dr and sublattice of Site1 and Site2: Sublattice[0] and Sublattice[1].
 */

class Distance{
public:
    int dSublattice;
    Vec<int> dCoordinate;
    
    Distance(){}
    Distance(const Site& i, const Site& j)
    {
        dCoordinate=j.Coordinate-i.Coordinate;
        dSublattice=GetSublatticeIndex(i.Sublattice, j.Sublattice);
    }
    
    Vec<real> GetVec();
    Distance Mirror();
    int GetSublattice(const int&);
    
    bool operator==(const Distance& v2)
    {
        if(dSublattice!=v2.dSublattice) return false;
        if(Dr!=v2.Dr) return false;
        return true;
    }
    
    bool operator!=(const Distance& v2)
    {
        if(*this==v2) return false;
        return true;
    }
};

Distance operator-(const Site& S1, const Site& S2);


/**
 *  class Lattice includes three set of vectors: 1) LatticeVec (unit cell lattice vector); 2)ReLatticeVec (reciprocal lattice vector for k); 3) SublatticeVec (vectors between different sublattices in the same unit cell).
 */
class Lattice{
public:
    Vec<real> LatticeVec[D];
    Vec<real> ReciprocalLatticeVec[D];
    Vec<real> SublatticeVec[NSublattice];
    
    Lattice() {Initialize();}
    void PlotLattice();
    
private:
    void Initialize();
};

extern Lattice lattice;

int Mirror(int, const int&);
int TestLattice();

#endif /* defined(__Fermion_Simulator__lattice__) */
