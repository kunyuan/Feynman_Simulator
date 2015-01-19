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
#include <string>
typedef double real;

#define M_RAN_INVM32 2.32830643653869628906e-010
#define M_RAN_INVM52 2.22044604925031308085e-016

/* Please refer to the following paper for more infomation on fast integer->float conversion
Jurgen A. Doornik. 2007. Conversion of high-period random numbers to floating point. ACM Trans. Model. Comput. Simul. 17, 1, Article 3 (January 2007). DOI=10.1145/1189756.1189759 http://doi.acm.org/10.1145/1189756.118975
*/
//double 32bits
#define RANDBL_32(iRan1) \
    ((int)(iRan1)*M_RAN_INVM32 + 0.5)

#define RANDBL_32_NO_ZERO(iRan1) \
    ((int)(iRan1)*M_RAN_INVM32 + (0.5 + M_RAN_INVM32 / 2))

//float number with 52bits
#define RANDBL_52_NO_ZERO(iRan1, iRan2) \
    ((int)(iRan1)*M_RAN_INVM32 + (0.5 + M_RAN_INVM52 / 2) + (int)((iRan2)&0x000FFFFF) * M_RAN_INVM52)

class RandomFactory {
    friend std::ostream& operator<<(std::ostream& os, RandomFactory& r);
    friend std::istream& operator>>(std::istream& is, RandomFactory& r);
    friend std::string ToString(const RandomFactory& RNG);

private:
    std::mt19937 _eng;

public:
    RandomFactory();
    RandomFactory(int);
    RandomFactory(const std::string& state);
    void Reset();
    void Reset(int);
    void Reset(const std::string& state);
    inline real urn()
    {
        return RANDBL_32(_eng());
    }

    //Generator integer random numbers in the closed interval [from, thru]
    inline int irn(int from, int thru)
    {
        return _eng() % (thru - from + 1) + from;
    }
    //    inline int nrn(real mean, real std)
    //    {
    //        static std::normal_distribution<> d{};
    //        return d(_eng, decltype(d)::param_type{mean, std});
    //    }
};

int TestRNG();
#endif
