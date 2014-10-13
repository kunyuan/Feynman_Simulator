//
//  rng.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/4/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "rng.h"

using namespace std;

RandomFactory::RandomFactory()
{
    random_device dev;
    _eng.seed(dev());
}

RandomFactory::RandomFactory(int seed)
{
    _eng.seed(seed);
}

void RandomFactory::Reset()
{
    random_device dev;
    _eng.seed(dev());
}

void RandomFactory::Reset(int seed)
{
    _eng.seed(seed);
}

real RandomFactory::urn()
{
    static uniform_real_distribution<real> d(0.0, 1.0);
    return d(_eng);
}

int RandomFactory::irn(int from, int thru)
{
    static std::uniform_int_distribution<> d{};
    return d(_eng, decltype(d)::param_type{from, thru});
}

int RandomFactory::irn1(int from, int thru)
{
    return _eng() % (thru - from + 1) + from;
}

// double pick_a_number(double from, double upto)
//{
//    static std::uniform_real_distribution<> d{};
//    using parm_t = decltype(d)::param_type;
//    return d(global_urng(), parm_t{from, upto});
//}

std::ostream &operator<<(std::ostream &os, RandomFactory &r)
{
    os << r._eng;
    return os;
}

std::istream &operator>>(std::istream &is, RandomFactory &r)
{
    is >> r._eng;
    return is;
}
