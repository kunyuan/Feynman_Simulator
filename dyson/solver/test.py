#!/usr/bin/env python
import unittest
import lu_fast as fast
import lu_slow as slow
import numpy as np

Shape=(32,32,100)
class TestWeightFFT(unittest.TestCase):
    def setUp(self):
        pass
    def test_slow(self):
        a=np.random.rand(*Shape)+1j*np.random.rand(*Shape)
        aa=a.copy()
        b=np.random.rand(*Shape)+1j*np.random.rand(*Shape)
        lu,piv=slow.lu_factor(a)
        x=slow.lu_solve(lu, piv, b)
        bb=np.einsum('jki,kli->jli', aa, x)
        self.assertTrue(np.allclose(b, bb, atol=1e-8))
    def test_fast(self):
        a=np.random.rand(*Shape)+1j*np.random.rand(*Shape)
        aa=a.copy()
        b=np.random.rand(*Shape)+1j*np.random.rand(*Shape)
        lu,piv=fast.lu_factor(a)
        x=fast.lu_solve(lu, piv, b)
        bb=np.einsum('jki,kli->jli', aa, x)
        self.assertTrue(np.allclose(b, bb, atol=1e-8))
    def test_det(self):
        a=np.random.rand(*Shape)+1j*np.random.rand(*Shape)
        #a=np.random.rand(*Shape)+1j*0.0
        aa=a.copy()
        aaa=a.copy()
        lu,piv=fast.lu_factor(a)
        #print np.diag(lu[:,:,0]), piv[:,0]
        x=fast.lu_det(lu, piv)
        lu,piv=slow.lu_factor(aa)
        #print np.diag(lu[:,:,0]), piv[:,0]
        xx=slow.lu_det(lu, piv)
        #print "det comparsion"
        #print x, xx, np.linalg.det(aaa[:,:,0])
        self.assertTrue(np.allclose(x, xx, atol=1e-6))

if __name__ == '__main__':
    unittest.main(verbosity=3)
