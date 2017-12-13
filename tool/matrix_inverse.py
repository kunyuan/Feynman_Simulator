import numpy as np
import timeit as t
B = np.random.rand(100000,3,3,2)
A=B[...,0]+1j*B[...,1]
#A = np.random.rand(100000,3,3)
def slow_inverse(A): 
    Ainv = np.zeros_like(A)
    for i in range(A.shape[0]):
        Ainv[i,:,:] = np.linalg.inv(A[i,:,:])
    return Ainv

def adjoint(A):
    """compute inverse without division by det; ...xv3xc3 input, or array of matrices assumed"""
    AI = np.empty_like(A)
    for i in xrange(3):
        AI[...,i,:] = np.cross(A[...,i-2,:], A[...,i-1,:])
    return AI

def inverse_transpose(A):
    """
    efficiently compute the inverse-transpose for stack of 3x3 matrices
    """
    I = adjoint(A)
    det = dot(I, A).mean(axis=-1)
    return I / det[...,None,None]

def inverse(A):
    """inverse of a stack of 3x3 matrices"""
    return np.swapaxes( inverse_transpose(A), -1,-2)
def dot(A, B):
    """dot arrays of vecs; contract over last indices"""
    return np.einsum('...i,...i->...', A, B)

def fast_inverse(A):
    identity = np.identity(A.shape[2], dtype=A.dtype)
    Ainv = np.zeros_like(A)

    for i in range(A.shape[0]):
        Ainv[i] = np.linalg.solve(A[i], identity)
    return Ainv

#def fast_inverse2(A):
    #identity = np.identity(A.shape[2], dtype=A.dtype)

    #return np.array([np.linalg.solve(x, identity) for x in A])

#from numpy.linalg import lapack_lite
#lapack_routine = lapack_lite.dgesv
## Looking one step deeper, we see that solve performs many sanity checks.  
## Stripping these, we have:
#def faster_inverse(A):
    #b = np.identity(A.shape[2], dtype=A.dtype)

    #n_eq = A.shape[1]
    #n_rhs = A.shape[2]
    #pivots = zeros(n_eq, np.intc)
    #identity  = np.eye(n_eq)
    #def lapack_inverse(a):
        #b = np.copy(identity)
        #pivots = zeros(n_eq, np.intc)
        #results = lapack_lite.dgesv(n_eq, n_rhs, a, n_eq, pivots, b, n_eq, 0)
        #if results['info'] > 0:
            #raise LinAlgError('Singular matrix')
        #return b

    #return array([lapack_inverse(a) for a in A])


x=100
print t.timeit("aI11 = slow_inverse(A)", setup="from __main__ import slow_inverse, A", number=1) 
I1=slow_inverse(A)
print np.dot(I1[x,:,:], A[x,:,:])
print t.timeit("aI11 = fast_inverse(A)", setup="from __main__ import fast_inverse, A", number=1) 
I2=fast_inverse(A)
print np.dot(I2[x,:,:], A[x,:,:])
print t.timeit("aI12 = inverse(A)", setup="from __main__ import inverse, A", number=1) 
I=inverse(A)
print np.dot(I[x,:,:], A[x,:,:])
