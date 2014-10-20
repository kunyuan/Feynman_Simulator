//
//  vector.cpp
//  Feynman_Simulator
//
//  Created by yuan on 10/14/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//
#include "vector.h"
#include <iostream>
#include <sstream>
#define SEP ' '

using namespace std;

template <typename T>
ostream &operator<<(ostream &os, Vec<T> &v)
{
    for (int i = 0; i < D - 1; i++)
        os << v[i] << SEP;
    os << v[D - 1];
    return os;
}

template ostream &operator<<(ostream &, Vec<int> &);
template ostream &operator<<(ostream &, Vec<real> &);

template <typename T>
istream &operator>>(istream &is, Vec<T> &v)
{
    for (int i = 0; i < D; i++)
        is >> v[i];
    return is;
}
template istream &operator>>(istream &, Vec<int> &);
template istream &operator>>(istream &, Vec<real> &);

template <typename T>
string Vec<T>::PrettyString()
{
    stringstream os;
    os << "(";
    for (int i = 0; i < D - 1; i++)
        os << _Arrary[i] << ",";
    os << _Arrary[D - 1] << ")";
    return os.str();
}

template <typename T>
int Vec<T>::ToIndex() const
{
    int Index = _Arrary[D - 1];
    for (int i = D - 2; i >= 0; i--) {
        Index = Index * L[i] + _Arrary[i];
    }
    return Index;
}

template class Vec<int>;
template class Vec<real>;
