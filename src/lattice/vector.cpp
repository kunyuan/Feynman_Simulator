//
//  vector.cpp
//  Feynman_Simulator
//
//  Created by yuan on 10/14/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//
#include "vector.h"
#include "utility.h"
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
        os << _Array[i] << ",";
    os << _Array[D - 1] << ")";
    return os.str();
}
template <typename T>
bool operator!=(const Vec<T>& v1, const Vec<T>&v2)
{
    if (v1 == v2)
        return false;
    return true;
}

template class Vec<int>;
template class Vec<real>;


bool operator==(const Vec<int>& v1, const Vec<int>&v2)
{
    for (int j = 0; j < D; j++)
        if (v1[j]!=v2[j])
            return false;
    return true;
}

bool operator==(const Vec<real>& v1, const Vec<real>&v2)
{
    for (int j = 0; j < D; j++)
        if (!Equal(v1[j],v2[j],eps0))
            return false;
    return true;
}