//
//  momentum.h
//  Feynman_Simulator
//
//  Created by yuan on 11/5/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__momentum__
#define __Feynman_Simulator__momentum__

#include <stdio.h>
#include <iostream>

const int MAX_K = 10000;

class Momentum{
  public:
    Momentum();// Default constructor
    Momentum(int k);
    Momentum(const Momentum &m); // Copy constructor

    int K; // member

    int index();
    int abs();
    
    Momentum &operator=(const Momentum &);

    // define the compound assignment operators first
    Momentum &operator+=(const Momentum &);
    Momentum &operator+=(int);
    Momentum &operator-=(const Momentum &);
    Momentum &operator-=(int);

    friend std::ostream &operator<<(std::ostream &, const Momentum &);
    friend std::istream &operator>>(std::istream &, Momentum &);
};

// Nonmember operators (to allow implicit conversion of the left operand)
Momentum operator+(const Momentum &, const Momentum &);
Momentum operator+(const Momentum &, int);
Momentum operator+(int, const Momentum &);

Momentum operator-(const Momentum &, const Momentum &);
Momentum operator-(const Momentum &, int);
Momentum operator-(int, const Momentum &);

Momentum operator*(int, const Momentum &);

bool operator==(const Momentum &, const Momentum &);
bool operator==(const Momentum &, int);
bool operator==(int, const Momentum &);

bool operator!=(const Momentum &, const Momentum &);
bool operator!=(const Momentum &, int);
bool operator!=(int, const Momentum &);

#endif /* defined(__Feynman_Simulator__momentum__) */
