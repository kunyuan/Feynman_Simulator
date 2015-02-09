import numpy as np
import scipy.linalg.lapack as lapack

def lu_factor(arr):
    """
       arr: complex type array with [M, N, N], it may be altered!!!
       return: (lu, piv): lu decomposition of array arr
    """
    if arr.shape[1] != arr.shape[2]:
        raise TypeError("the input array must have shape [:, N, N]")
    piv=np.zeros((arr.shape[0],arr.shape[1]))
    arr=np.ascontiguousarray(arr)
    for index in range(arr.shape[0]):
        arr[index,:,:],piv[index,:],info=lapack.zgetrf(arr[index,:,:],overwrite_a=True)
        if info < 0:
            raise ValueError('illegal value in %d-th argument of internal getrf' % -info)
        elif info > 0:
            raise ValueError("Diagonal number %d is exactly zero. Singular matrix." % -info)
    return arr, piv

def lu_solve(lu, piv , b):
    """
       lu, piv: lu decomposition of array arr
       b: shape [M, N, N], it will be altered!!!
       return: the solve of aX=b, an array with shape [M,N,N]
    """
    if lu.shape[1] != lu.shape[2]:
        raise TypeError("the lu array must have shape [:, N, N]")
    if lu.shape[0] != piv.shape[0]:
        raise TypeError("require lu.shape[0]==piv.shape[0]")
    if lu.shape[1] != piv.shape[1]:
        raise TypeError("require lu.shape[1]==piv.shape[1]")
    if lu.shape[0] != b.shape[0] or lu.shape[1]!=b.shape[1] or lu.shape[2]!=b.shape[2]:
        raise TypeError("require lu.shape==b.shape")
    """solv ax=b, self.Data will be from lu of a to x"""
    lu=np.ascontiguousarray(lu)
    piv=np.ascontiguousarray(piv)
    b=np.ascontiguousarray(b)
    for index in range(lu.shape[0]):
        b[index,:,:],info = lapack.zgetrs(lu[index,:,:], piv[index,:], b[index,:,:], overwrite_b=False)
        if info is not 0:
            raise ValueError('illegal value in %d-th argument of internal gesv|posv'% -info)
    return b

def lu_det(lu, piv):
    if lu.shape[1] != lu.shape[2]:
        raise TypeError("the lu array must have shape [:, N, N]")
    if lu.shape[0] != piv.shape[0]:
        raise TypeError("require lu.shape[0]==piv.shape[0]")
    if lu.shape[1] != piv.shape[1]:
        raise TypeError("require lu.shape[1]==piv.shape[1]")
    Det=np.ones(lu.shape[0], dtype=lu.dtype)
    diagrange=range(lu.shape[1])
    for index in range(lu.shape[0]):
        for i in diagrange:
            if piv[index, i] is not i:
                Det[index] *=-lu[index,i,i]
            else:
                Det[index] *=lu[index,i,i]
    return Det

