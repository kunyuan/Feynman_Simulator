//---------------------------------------------------------------------------
// �Nikolai V. Shokhirev, 2004-2008  <nikolai@shokhirev.com>
// http://www.shokhirev.com/nikolai.html
// Reduced demo version
//---------------------------------------------------------------------------

#ifndef MathUtils_H
#define MathUtils_H

#include <math.h>
#include <string>
#include <sstream>

#ifndef real
#define real double
#endif

//---------------------------------------------------------------------------

const real eps0 = 1.0e-9;

// Macheps + 1.0 > 1.0
const real Macheps = 2.22044604925031E-16; // for double

const real MaxReal = 1.0e30;
const real MinReal = -1.0e30;

// FORTRAN abs
real abs(real x); // { return ( (x >= 0.0) ? x : -x); }

// FORTRAN iabs
int iabs(int x); // { return ( (x >= 0.0) ? x : -x); }
// more functions ...

//   if ( abs(x)> eps) return x; else return 0.0;
real Zero(real x, real eps);

// float equal
bool Equal(real x1, real x2, real eps = eps0);

// initial an array with some value
template <typename T>
static void InitialArry(T *array, int num, T value)
{
    for (int i = 0; i < num; i++) {
        *(array + i) = value;
    }
}

template <typename T>
std::string ToString(const T &value)
{
    std::ostringstream oss;
    oss << value;
    return oss.str();
}

template <typename T>
void AssignFromTo(T *source, T *target, int size)
{
    for (int i = 0; i < size; i++)
        source[i] = target[i];
}

#endif
