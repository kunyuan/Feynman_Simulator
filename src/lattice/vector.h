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
#include <iostream>
#include "convention.h"
using namespace std;

template <typename T>
class Vec {
  private:
    T _Array[D];

  public:
    Vec()
    {
    }
    
    Vec(T value)
    {
        for(int j=0; j<D; j++)
            _Array[j] = value;
    }


    T &operator[](int index)
    {
        return _Array[index];
    }

    const T &operator[](int index) const
    {
        return _Array[index];
    }

    Vec operator*(int i) const
    {
        Vec v;
        for (int j = 0; j < D; j++)
            v[j] = _Array[j] * i;
        return v;
    }

    Vec operator*(const real &i) const
    {
        Vec v;
        for (int j = 0; j < D; j++)
            v[j] = _Array[j] * i;
        return v;
    }

    Vec operator+(const Vec &v2) const
    {
        Vec v;
        for (int j = 0; j < D; j++)
            v[j] = _Array[j] + v2._Array[j];
        return v;
    }

    Vec operator-(const Vec &v2) const
    {
        Vec v;
        for (int j = 0; j < D; j++)
            v[j] = _Array[j] - v2._Array[j];
        return v;
    }

    Vec &operator+=(const Vec &v2)
    {
        for (int j = 0; j < D; j++)
            _Array[j] += v2._Array[j];
        return *this;
    }

    Vec &operator-=(const Vec &v2)
    {
        for (int j = 0; j < D; j++)
            _Array[j] -= v2._Array[j];
        return *this;
    }

    string PrettyString();


    template <typename TT>
    friend std::ostream &operator<<(std::ostream &os, Vec<TT> &);
    template <typename TT>
    friend std::istream &operator>>(std::istream &is, Vec<TT> &);
    
    friend bool operator==(const Vec<int>&, const Vec<int>&);
    friend bool operator==(const Vec<real>&, const Vec<real>&);
    template <typename TT>
    friend bool operator!=(const Vec<TT>&, const Vec<TT>&);
};

#endif /* defined(__Feynman_Simulator__vector__) */
