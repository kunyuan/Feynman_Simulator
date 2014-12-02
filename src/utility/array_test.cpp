#include <math.h>

// Turning off argument checking with -DNDEBUG improves optimization.

#include "array.h"
using namespace Array;
using namespace std;

int n = 2, m = 3, p = 4;

typedef array1<double>::opt vector;

using std::cout;

void f(double *x)
{
    cout << x[0] << endl;
    return;
}

template <class T>
void g(T x)
{
    cout << x << endl;
    return;
}

template <class T>
double h(typename array1<T>::opt &x)
{
    return x[0];
}

int TestArray()
{
    array3<double> A(n, m, p);
    double sum = 0.0;

    // Sequential access:
    int size = A.Size();
    for (int i = 0; i < size; i++)
        A(i) = i;

    // Random access:
    for (int i = 0; i < n; i++) {

        // The following statements are equivalent, but the first one optimizes better.
        array2<double> Ai = A[i];
        //Ai is just an array with temporary state, the ownwer of data is still A
        //        array2<double> Ai;
        //        Ai = A[i];
        //    array2<double> Ai(m,p); Ai=A[i]; // This does an extra memory copy.

        for (int j = 0; j < m; j++) {
            //      array1<double> Aij=Ai[j];
            // For 1D arrays: many compilers optimize array1<>::opt better than array1<>.
            vector Aij = Ai[j];

            for (int k = 0; k < p; k++) {
                // The following statements are equivalent, but the first one optimizes better.
                sum = sum + Aij[k];
                //	sum=sum+A(i,j,k); // This does extra index multiplication.
            }
        }
    }

    cout << sum << endl;

    f(A);
    g(A);

    array2<int> temp;
    temp.Allocate(1, 2);
    //the technique may be used to save one destruction of the array data

    return 0;
}
