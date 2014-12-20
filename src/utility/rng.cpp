//
//  rng.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/4/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "rng.h"
#include <sstream>

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

RandomFactory::RandomFactory(const std::string& state)
{
    Reset(state);
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

void RandomFactory::Reset(const std::string& state)
{
    std::istringstream iss(state);
    iss >> _eng;
}

// double pick_a_number(double from, double upto)
//{
//    static std::uniform_real_distribution<> d{};
//    using parm_t = decltype(d)::param_type;
//    return d(global_urng(), parm_t{from, upto});
//}

std::ostream& operator<<(std::ostream& os, RandomFactory& r)
{
    os << r._eng;
    return os;
}

std::istream& operator>>(std::istream& is, RandomFactory& r)
{
    is >> r._eng;
    return is;
}

std::string ToString(const RandomFactory& rng)
{
    std::ostringstream oss;
    oss << rng._eng;
    return oss.str();
}
RandomFactory RNG;
