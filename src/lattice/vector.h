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
    T _Arrary[D];

  public:
    Vec()
    {
    }

    Vec(int index)
    {
        for (int i = 0; i < D - 1; i++) {
            _Arrary[i] = index % L[i];
            index = index / L[i];
        }
        _Arrary[D - 1] = index;
    }

    T &operator[](int index)
    {
        return _Arrary[index];
    }

    const T &operator[](int index) const
    {
        return _Arrary[index];
    }

    Vec operator*(int i) const
    {
        Vec v;
        for (int j = 0; j < D; j++)
            v[j] = _Arrary[j] * i;
        return v;
    }

    Vec operator*(const real &i) const
    {
        Vec v;
        for (int j = 0; j < D; j++)
            v[j] = _Arrary[j] * i;
        return v;
    }

    Vec operator+(const Vec &v2) const
    {
        Vec v;
        for (int j = 0; j < D; j++)
            v[j] = _Arrary[j] + v2._Arrary[j];
        return v;
    }

    Vec operator-(const Vec &v2) const
    {
        Vec v;
        for (int j = 0; j < D; j++)
            v[j] = _Arrary[j] - v2._Arrary[j];
        return v;
    }

    Vec &operator+=(const Vec &v2)
    {
        for (int j = 0; j < D; j++)
            _Arrary[j] += v2._Arrary[j];
        return *this;
    }

    Vec &operator-=(const Vec &v2)
    {
        for (int j = 0; j < D; j++)
            _Arrary[j] -= v2._Arrary[j];
        return *this;
    }

    bool operator==(const Vec &v2) const
    {
        for (int j = 0; j < D; j++)
            if (_Arrary[j] != v2._Arrary[j])
                return false;
        return true;
    }

    bool operator!=(const Vec &v2) const
    {
        if (*this == v2)
            return false;
        return true;
    }

    string PrettyString();

    int ToIndex() const;

    template <typename TT>
    friend std::ostream &operator<<(std::ostream &os, Vec<TT> &);
    template <typename TT>
    friend std::istream &operator>>(std::istream &is, Vec<TT> &);
};

int TestVector();
#endif /* defined(__Feynman_Simulator__vector__) */
