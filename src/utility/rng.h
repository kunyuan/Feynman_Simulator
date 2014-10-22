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

class RandomFactory {
    friend std::ostream &operator<<(std::ostream &os, RandomFactory &r);
    friend std::istream &operator>>(std::istream &is, RandomFactory &r);

  private:
    std::mt19937 _eng;

  public:
    RandomFactory();
    RandomFactory(int);
    void Reset();
    void Reset(int);
    inline real urn()
    {
        static std::uniform_real_distribution<real> d(0.0, 1.0);
        return d(_eng);
    }
    //Generator integer random numbers in the closed interval [from, thru]
    inline int irn(int from, int thru)
    {
        static std::uniform_int_distribution<> d{};
        return d(_eng, decltype(d)::param_type{from, thru});
    }

    //Generator integer random numbers in the closed interval [from, thru]
    inline int irn1(int from, int thru)
    {
        return _eng() % (thru - from + 1) + from;
    }
    //    inline int nrn(real mean, real std)
    //    {
    //        static std::normal_distribution<> d{};
    //        return d(_eng, decltype(d)::param_type{mean, std});
    //    }
};

extern RandomFactory RNG;

int TestRNG();
#endif
