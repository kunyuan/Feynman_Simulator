//
//  vector.h
//  Feynman_Simulator
//
//  Created by yuan on 10/14/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__vector__
#define __Feynman_Simulator__vector__

#include <stdio.h>
#include "convention.h"

template <typename T>
class Vec{
private:
    T _Arrary[D];
public:
    Vec()
    {
        for(int i=0; i<D; i++)
            _Arrary[i]=0;
    }
    
    Vec(T t)
    {
        for(int i=0; i<D; i++)
            _Arrary[i]=t;
    }
    
    T& operator[](const int& index)
    {
        return _Arrary[index];
    }
    
    Vec operator*(const int& i) const
    {
        Vec v;
        for(int j=0; j<D; j++)
            v[j]=_Arrary[j]*i;
        return v;
    }

    Vec operator*(const real& i) const
    {
        Vec v;
        for(int j=0; j<D; j++)
            v[j]=_Arrary[j]*i;
        return v;
    }
    
    Vec operator+(const Vec& v2) const
    {
        Vec v;
        for(int j=0; j<D; j++)
            v[j]=_Arrary[j]+v2._Arrary[j];
        return v;
    }
    
    Vec& operator+=(const Vec& v2)
    {
        for(int j=0; j<D; j++)
            _Arrary[j]+=v2._Arrary[j];
        return *this;
    }
    
    bool operator==(const Vec& v2) const
    {
        for(int j=0; j<D; j++)
            if(_Arrary[j]!=v2._Arrary[j]) return false;
        return true;
    }
};

int TestVector();
#endif /* defined(__Feynman_Simulator__vector__) */
