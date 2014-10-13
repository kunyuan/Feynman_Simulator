//---------------------------------------------------------------------------
// �Nikolai V. Shokhirev, 2004-2008  <nikolai@shokhirev.com>
// http://www.shokhirev.com/nikolai.html
// Reduced demo version
//---------------------------------------------------------------------------

#include "utility.h"

using namespace std;

//---------------------------------------------------------------------------

// FORTRAN abs
real abs(real x)
{
    return ((x >= 0.0) ? x : -x);
}

// FORTRAN iabs
int iabs(int x)
{
    return ((x >= 0.0) ? x : -x);
}
// more functions ...

//   if ( abs(x)> eps) return x; else return 0.0;
real Zero(real x, real eps)
{
    return ((abs(x) > eps) ? x : 0.0);
}

// float equal
bool Equal(real x1, real x2, real eps)
{
    return (fabs(x1 - x2) < eps);
};
