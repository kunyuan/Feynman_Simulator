//
//  vector.h
//  Feynman_Simulator
//
//  Created by yuan on 10/14/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__vector__
#define __Feynman_Simulator__vector__

#include <iosfwd>
#include <sstream>
#include <vector>
#include <initializer_list>
#include "../utility/convention.h"
using namespace std;

template <typename T>
class Vec {
private:
    T _Array[D];

public:
    Vec()
    {
    }

    Vec(std::initializer_list<T> list)
    {
        int i = 0;
        for (auto p = list.begin(); p < list.end() && i < D; ++p) {
            _Array[i] = *p;
            i++;
        }
    }

    Vec(T value)
    {
        for (int j = 0; j < D; j++)
            _Array[j] = value;
    }

    Vec(T* value)
    {
        for (int j = 0; j < D; j++)
            _Array[j] = value[j];
    }

    void CopyToArray(T* target) const
    {
        for (int i = 0; i < D; i++)
            target[i] = _Array[i];
    }

    const T* begin() const
    {
        return _Array;
    }
    const T* end() const
    {
        return _Array + D;
    }
    uint size() const
    {
        return D;
    }

    T& operator[](int index)
    {
        return _Array[index];
    }

    const T& operator[](int index) const
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

    Vec operator*(const real& i) const
    {
        Vec v;
        for (int j = 0; j < D; j++)
            v[j] = _Array[j] * i;
        return v;
    }

    Vec operator+(const Vec& v2) const
    {
        Vec v;
        for (int j = 0; j < D; j++)
            v[j] = _Array[j] + v2._Array[j];
        return v;
    }

    Vec operator-(const Vec& v2) const
    {
        Vec v;
        for (int j = 0; j < D; j++)
            v[j] = _Array[j] - v2._Array[j];
        return v;
    }

    Vec& operator+=(const Vec& v2)
    {
        for (int j = 0; j < D; j++)
            _Array[j] += v2._Array[j];
        return *this;
    }

    Vec& operator-=(const Vec& v2)
    {
        for (int j = 0; j < D; j++)
            _Array[j] -= v2._Array[j];
        return *this;
    }

    string PrettyString();

    template <typename TT>
    friend std::ostream& operator<<(std::ostream& os, const Vec<TT>&);
    template <typename TT>
    friend std::istream& operator>>(std::istream& is, Vec<TT>&);

    friend bool operator==(const Vec<int>&, const Vec<int>&);
    friend bool operator==(const Vec<real>&, const Vec<real>&);
};
bool operator!=(const Vec<int>& v1, const Vec<int>& v2);
bool operator!=(const Vec<real>& v1, const Vec<real>& v2);

template <typename T>
std::string ToString(Vec<T> value)
{
    std::ostringstream oss;
    oss << value;
    return oss.str();
}

#endif /* defined(__Feynman_Simulator__vector__) */
