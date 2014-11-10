//
//  momentum.cpp
//  Feynman_Simulator
//
//  Created by yuan on 11/5/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "momentum.h"
#include "utility/rng.h"

Momentum::Momentum():K(0)// Default constructor
{
}

Momentum::Momentum(int k):K(k)
{
}

Momentum::Momentum(const Momentum &m):K(m.K) // Copy constructor
{
}

Momentum& Momentum::operator=(int k)
{
    if(K!=k)  K = k;
    return (*this);
}

Momentum& Momentum::operator=(const Momentum& m)
{
    if(K!=m.K)  K = m.K;
    return (*this);
}


Momentum& Momentum::operator+=(int k)
{
    K += k;
    if(K>=MAX_K) K -= MAX_K;
    return (*this);
}

Momentum& Momentum::operator+=(const Momentum& m)
{
    K += m.K;
    if(K>=MAX_K) K -= MAX_K;
    return (*this);
}

Momentum& Momentum::operator-=(int k)
{
    K -= k;
    if(K<-MAX_K) K += MAX_K;
    return (*this);
}

Momentum& Momentum::operator-=(const Momentum& m)
{
    K -= m.K;
    if(K<-MAX_K) K += MAX_K;
    return (*this);
}


Momentum operator+(const Momentum & m1, int k)
{
    return (Momentum)(m1.K+k);
}

Momentum operator+(const Momentum & m1, const Momentum & m2)
{
    return m1+m2.K;
}

Momentum operator+(int k, const Momentum & m1)
{
    return m1+k;
}

Momentum operator-(const Momentum & m1, int k)
{
    return (Momentum)(m1.K-k);
}

Momentum operator-(const Momentum & m1, const Momentum & m2)
{
    return m1-m2.K;
}

Momentum operator-(int k, const Momentum & m1)
{
    return (Momentum)(k-m1.K);
}

Momentum operator*(int k, const Momentum & m1)
{
    return (Momentum)(m1.K*k);
}

bool operator==(const Momentum & m1, const Momentum & m2)
{
    return (m1.K==m2.K);
}

bool operator==(const Momentum & m1, int k)
{
    return (m1.K==k);
}

bool operator==(int k, const Momentum & m1)
{
    return (m1.K==k);
}

bool operator!=(const Momentum & m1, const Momentum & m2)
{
    return (m1.K!=m2.K);
}

bool operator!=(const Momentum & m1, int k)
{
    return (m1.K!=k);
}

bool operator!=(int k, const Momentum & m1)
{
    return (m1.K!=k);
}

std::ostream & operator<<(std::ostream & os, const Momentum & m)
{
    os << m.K;
    return os;
}

std::istream & operator>>(std::istream & is, Momentum & m)
{
    is >> m.K;
    return is;
}


