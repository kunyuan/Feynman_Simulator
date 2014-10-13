//
//  random.h
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/3/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef Fermion_Simulator_random_h
#define Fermion_Simulator_random_h

#include <random>
#ifndef real
#define real double
#endif

using namespace std;

class RandomFactory {
    friend std::ostream &operator<<(std::ostream &os, RandomFactory &r);
    friend std::istream &operator>>(std::istream &is, RandomFactory &r);

  private:
    mt19937 _eng;

  public:
    RandomFactory();

    RandomFactory(int);
    void Reset();
    void Reset(int);
    real urn();
    int irn(int, int);
    int irn1(int, int);
};

int TestRNG();
#endif
