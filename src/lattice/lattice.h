//
//  lattice.h
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/6/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Fermion_Simulator__lattice__
#define __Fermion_Simulator__lattice__

const int NSublattice=2;
const int L[3]={4,4,1};
const int Vol=L[0]*L[1]*L[2];
const int D=2;

int a();

template <typename T>
class Vec{
private:
    T _Arrary[D];
public:
    Vec(T t)
    {
        for(int i=0; i<D; i++)
            _Arrary[i]=t;
    }
    
    T& operator[](int index)
    {
        return _Arrary[index];
    }
};

//class Lattice{
//};

//class Site{
//public:
//    Site();
//    bool Valid() const;
//    int Sublattice;
//};

//int Mirror(const int&, const int&);

int TestLattice();
#endif /* defined(__Fermion_Simulator__lattice__) */
