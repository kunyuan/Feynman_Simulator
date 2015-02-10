import numpy as np
import scipy.linalg.lapack as lapack

def lu_factor(arr):
    """
       arr: complex type array with [N, N, M], it may be altered!!!
       return: (lu, piv): lu decomposition of array arr
    """
    if arr.shape[0] != arr.shape[1]:
        raise TypeError("the input array must have shape [N, N, :]")
    piv=np.zeros((arr.shape[0], arr.shape[2]), dtype=np.int, order=['F'])
    #arr=np.ascontiguousarray(arr)
    for index in range(arr.shape[2]):
        arr[:,:,index],piv[:,index],info=lapack.zgetrf(arr[:,:,index],overwrite_a=True)
        if info < 0:
            raise ValueError('illegal value in %d-th argument of internal getrf' % -info)
        elif info > 0:
            raise ValueError("Diagonal number %d is exactly zero. Singular matrix." % -info)
    return arr, piv

def lu_solve(lu, piv , b):
    """
       lu, piv: lu decomposition of array arr
       b: shape [M, N, N], it will not be altered
       return: the solve of aX=b, an array with shape [M,N,N]
    """
    if lu.shape[0] != lu.shape[1]:
        raise TypeError("the lu array must have shape [N, N, :]")
    if lu.shape[2] != piv.shape[1]:
        raise TypeError("require lu.shape[2]==piv.shape[1]")
    if lu.shape[1] != piv.shape[0]:
        raise TypeError("require lu.shape[1]==piv.shape[0]")
    if lu.shape[2] != b.shape[2] or lu.shape[1]!=b.shape[1] or lu.shape[0]!=b.shape[0]:
        raise TypeError("require lu.shape==b.shape")
    """solv ax=b, self.Data will be from lu of a to x"""
    #lu=np.ascontiguousarray(lu)
    #piv=np.ascontiguousarray(piv)
    bb=np.ascontiguousarray(b)
    if np.may_share_memory(b,bb):
        bb=b.copy()
    for index in range(lu.shape[2]):
        bb[:,:,index],info = lapack.zgetrs(lu[:,:,index], piv[:,index], bb[:,:,index], overwrite_b=False)
        if info is not 0:
            raise ValueError('illegal value in %d-th argument of internal gesv|posv'% -info)
    return bb

def lu_det(lu, piv):
    if lu.shape[0] != lu.shape[1]:
        raise TypeError("the lu array must have shape [N, N,:]")
    if lu.shape[2] != piv.shape[1]:
        raise TypeError("require lu.shape[2]==piv.shape[1]")
    if lu.shape[1] != piv.shape[0]:
        raise TypeError("require lu.shape[1]==piv.shape[0]")
    Det=np.ones(lu.shape[2], dtype=lu.dtype)
    diagrange=range(lu.shape[1])
    for index in range(lu.shape[2]):
        for i in diagrange:
            if piv[i,index]!=i:
                Det[index] *=-lu[i,i,index]
            else:
                Det[index] *=lu[i,i,index]
    return Det

