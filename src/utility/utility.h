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
#include "convention.h"
#include "abort.h"

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

//float iszero
bool Zero(real x, real eps = eps0);

// float equal
bool Equal(real x1, real x2, real eps = eps0);
bool Equal(uint x1, uint x2, real eps = eps0);
template <typename T>
bool Equal(const T* x1, const T* x2, uint num, real eps = eps0)
{
    for (uint i = 0; i < num; i++)
        if (!Equal(x1[i], x2[i], eps0))
            return false;
    return true;
}

template <typename T>
std::string ToString(const T& value)
{
    std::ostringstream oss;
    oss << value;
    return oss.str();
}

/* Convert double to string with specified number of places after the decimal
   and left padding. */
std::string ToString(const double x, const int width, const int decDigits = 5);
/*! Center-aligns string within a field of width w. Pads with blank spaces
    to enforce alignment. */
std::string Center(const std::string s, const int w);
/* Right-aligns string within a field of width w. Pads with blank spaces
 to enforce alignment. */
std::string Right(const std::string s, const int w);
/*! Left-aligns string within a field of width w. Pads with blank spaces
 to enforce alignment. */
std::string Left(const std::string s, const int w);

template <typename T>
void AssignFromTo(T* source, T* target, int size)
{
    for (int i = 0; i < size; i++)
        target[i] = source[i];
}

template <typename T>
void InitialArray(T* target, T t, const int& size)
{
    for (int i = 0; i < size; i++)
        target[i] = t;
}

bool CleanFile(const std::string& FileName);

bool DoesFileExist(const std::string& FileName);

#define CHECKNULL(source)                     \
    {                                         \
        if ((source) == nullptr)              \
            ABORT(#source << " is nullptr!"); \
    }

#endif

