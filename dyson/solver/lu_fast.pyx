import numpy as np
cimport numpy as np
import cython
from cython.parallel import prange

ctypedef np.double_t cDOUBLE
ctypedef np.complex128_t cComplex
DOUBLE = np.float64
Complex = np.complex128

#detail comments in the end
cdef extern void zgetrf_(int* M, int* N, cComplex* a, int* lda, int* ipiv, int* info);
#*  ZGETRF computes an LU factorization of a general M-by-N matrix A
#*  using partial pivoting with row interchanges.
cdef extern void zgetrs_(char *TRANS, int *N, int *NRHS, cComplex *A, int *LDA, int *IPIV, cComplex *B, int *LDB, int *INFO );
#*  ZGETRS solves a system of linear equations
#*     A * X = B,  A**T * X = B,  or  A**H * X = B
#*  with a general N-by-N matrix A using the LU factorization computed
#*  by ZGETRF.

def lu_factor(np.ndarray[cComplex, ndim=3] arr):
    """
       arr: complex type array with [N, N, M], it may be altered!!!
       return: (lu, piv): lu decomposition of array arr in fortran order!!!
    """
    if arr.shape[0] != arr.shape[1]:
        raise TypeError("the input array must have shape [N, N, :]")
    cdef int size= arr.shape[2];
    cdef int N = arr.shape[0];
    cdef int LDA = N;
    cdef int info;
    cdef int i;
    cdef np.ndarray[int, ndim=2] piv;
    piv=np.empty((N, size), dtype=np.intc, order=['F'])
    #make arr fortran contiguous!
    arr=np.require(arr, requirements=['F'])
    #for i in prange(size, nogil=True):
    for i in range(size):
        zgetrf_(&N, &N, &arr[0,0,i], &LDA, &piv[0,i], &info);
        if info < 0:
            raise ValueError('illegal value in %d-th argument of internal getrf' % -info)
        elif info > 0:
            raise ValueError("Diagonal number %d is exactly zero. Singular matrix." % -info)
    #piv stores pivot indices with fortran convention, meaning indices start from 1
    return arr,piv-1

def lu_solve(np.ndarray[cComplex, ndim=3] lu, np.ndarray[int, ndim=2] piv, np.ndarray[cComplex, ndim=3] b):
    """
       lu, piv: lu decomposition of array arr
       b: shape [N, N, M], it will not be altered
       return: the solve of aX=b, an array with shape [N,N,M], but in fortran order!!!
    """
    if lu.shape[0] != lu.shape[1]:
        raise TypeError("the lu array must have shape [N, N, :]")
    if lu.shape[2] != piv.shape[1]:
        raise TypeError("require lu.shape[2]==piv.shape[1]")
    if lu.shape[1] != piv.shape[0]:
        raise TypeError("require lu.shape[1]==piv.shape[0]")
    if lu.shape[2] != b.shape[2] or lu.shape[1]!=b.shape[1] or lu.shape[0]!=b.shape[0]:
        raise TypeError("require lu.shape==b.shape")
    cdef char trans = 'N';
    cdef int size= lu.shape[2];
    cdef int N = lu.shape[0];
    cdef int nrhs = N;
    cdef int LDA = N;
    cdef int LDB = N;
    cdef int info;
    cdef np.ndarray[cComplex, ndim=3] bb
    lu=np.require(lu, requirements=['F'])
    piv=np.require(piv, requirements=['F'])+1
    bb=np.require(b, requirements=['F'])
    if np.may_share_memory(b,bb):
        bb=b.copy()
    #for i in prange(size, nogil=True):
    for i in range(size):
        zgetrs_(&trans, &N, &nrhs, &lu[0,0,i], &LDA, &piv[0,i], &bb[0,0,i], &LDB, &info);
        if info is not 0:
            raise ValueError('illegal value in %d-th argument of internal gesv|posv'% -info)
    return bb

#@cython.wraparound(False)
#@cython.boundscheck(False)
#@cython.nonecheck(False)
def lu_det(np.ndarray[cComplex, ndim=3] lu, np.ndarray[int, ndim=2] piv):
    if lu.shape[0] != lu.shape[1]:
        raise TypeError("the lu array must have shape [N, N, M]")
    if lu.shape[2] != piv.shape[1]:
        raise TypeError("require lu.shape[2]==piv.shape[1]")
    if lu.shape[1] != piv.shape[0]:
        raise TypeError("require lu.shape[1]==piv.shape[0]")
    cdef np.ndarray[cComplex, ndim=1] det
    cdef int size= lu.shape[2];
    cdef int N = lu.shape[0];
    #print lu.shape[2], lu.shape[0], N
    det=np.ones(size, dtype=Complex)
    for i in range(size):
        for j in range(N):
            if piv[j,i]!=j:
                det[i]=-det[i]*lu[j,j,i]
            else:
                det[i]=det[i]*lu[j,j,i]
            #print j,piv[j,i], lu[j,j,i], det[i]
    return det

#*     SUBROUTINE ZGETRF( M, N, A, LDA, IPIV, INFO )
#*
#*  -- LAPACK routine (version 3.1) --
#*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
#*     November 2006
#*
#*     .. Scalar Arguments ..
      #INTEGER            INFO, LDA, M, N
#*     ..
#*     .. Array Arguments ..
      #INTEGER            IPIV( * )
      #COMPLEX*16         A( LDA, * )
#*     ..
#*
#*  Purpose
#*  =======
#*
#*  ZGETRF computes an LU factorization of a general M-by-N matrix A
#*  using partial pivoting with row interchanges.
#*
#*  The factorization has the form
#*     A = P * L * U
#*  where P is a permutation matrix, L is lower triangular with unit
#*  diagonal elements (lower trapezoidal if m > n), and U is upper
#*  triangular (upper trapezoidal if m < n).
#*
#*  This is the right-looking Level 3 BLAS version of the algorithm.
#*
#*  Arguments
#*  =========
#*
#*  M       (input) INTEGER
#*          The number of rows of the matrix A.  M >= 0.
#*
#*  N       (input) INTEGER
#*          The number of columns of the matrix A.  N >= 0.
#*
#*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
#*          On entry, the M-by-N matrix to be factored.
#*          On exit, the factors L and U from the factorization
#*          A = P*L*U; the unit diagonal elements of L are not stored.
#*
#*  LDA     (input) INTEGER
#*          The leading dimension of the array A.  LDA >= max(1,M).
#*
#*  IPIV    (output) INTEGER array, dimension (min(M,N))
#*          The pivot indices; for 1 <= i <= min(M,N), row i of the
#*          matrix was interchanged with row IPIV(i).
#*
#*  INFO    (output) INTEGER
#*          = 0:  successful exit
#*          < 0:  if INFO = -i, the i-th argument had an illegal value
#*          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
#*                has been completed, but the factor U is exactly
#*                singular, and division by zero will occur if it is used
#*                to solve a system of equations.
#*
#*  =====================================================================
#*
      #SUBROUTINE ZGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
#*
#*  -- LAPACK routine (version 3.1) --
#*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
#*     November 2006
#*
#*     .. Scalar Arguments ..
      #CHARACTER          TRANS
      #INTEGER            INFO, LDA, LDB, N, NRHS
#*     ..
#*     .. Array Arguments ..
      #INTEGER            IPIV( * )
      #COMPLEX*16         A( LDA, * ), B( LDB, * )
#*     ..
#*
#*  Purpose
#*  =======
#*
#*  ZGETRS solves a system of linear equations
#*     A * X = B,  A**T * X = B,  or  A**H * X = B
#*  with a general N-by-N matrix A using the LU factorization computed
#*  by ZGETRF.
#*
#*  Arguments
#*  =========
#*
#*  TRANS   (input) CHARACTER*1
#*          Specifies the form of the system of equations:
#*          = 'N':  A * X = B     (No transpose)
#*          = 'T':  A**T * X = B  (Transpose)
#*          = 'C':  A**H * X = B  (Conjugate transpose)
#*
#*  N       (input) INTEGER
#*          The order of the matrix A.  N >= 0.
#*
#*  NRHS    (input) INTEGER
#*          The number of right hand sides, i.e., the number of columns
#*          of the matrix B.  NRHS >= 0.
#*
#*  A       (input) COMPLEX*16 array, dimension (LDA,N)
#*          The factors L and U from the factorization A = P*L*U
#*          as computed by ZGETRF.
#*
#*  LDA     (input) INTEGER
#*          The leading dimension of the array A.  LDA >= max(1,N).
#*
#*  IPIV    (input) INTEGER array, dimension (N)
#*          The pivot indices from ZGETRF; for 1<=i<=N, row i of the
#*          matrix was interchanged with row IPIV(i).
#*
#*  B       (input/output) COMPLEX*16 array, dimension (LDB,NRHS)
#*          On entry, the right hand side matrix B.
#*          On exit, the solution matrix X.
#*
#*  LDB     (input) INTEGER
#*          The leading dimension of the array B.  LDB >= max(1,N).
#*
#*  INFO    (output) INTEGER
#*          = 0:  successful exit
#*          < 0:  if INFO = -i, the i-th argument had an illegal value
#*
#*  =====================================================================
