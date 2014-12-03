#!/usr/bin/env python
import math
import numpy as np
import sys
import os.path
import unittest
from logger import *

SPIN,SPIN2,SPIN3=2,4,8
IN,OUT=0,1
DOWN,UP=0,1
SP,SUB,VOL,TAU=0,1,2,3
NSublattice=2

SHAPE={SUB: NSublattice**2, VOL: None, TAU : None}
def CheckShape(shape_):
    for e in SHAPE:
        if e>=len(shape_):
            pass
        elif SHAPE[e] is None:
            SHAPE[e]=shape_[e]
        else:
            Assert(SHAPE[e]==shape_[e], "Shape {0} does not match {1} at the {2}rd element!".format(SHAPE, shape_,e))

class IndexMap:
    def __init__(self,Beta, L_):
        self.__Beta=Beta
        self.__dBeta=Beta/SHAPE[TAU]
        self.__dBetaInverse=1.0/self.__dBeta
        self.__L=L_

    def TauIndex(self, In, Out):
        tau=Out-In
        Index=math.floor(tau*self.__dBetaInverse)
        return Index if tau>=0 else Index+SHAPE[TAU]

    def IndexToTau(self, Index):
        return (Index+0.5)*self.__dBeta

    def SublatIndex(self, In, Out):
        return In*NSublattice+Out
    def CoordiIndex(self,In, Out):
        #Out[0]*L1*L2+Out[1]*L2+Out[2] with In=(0,0,0)
        Index=Out[0]-In[0]
        for i in range(1,len(Vec)):
            Index=Index*self.__L[i]+(Out[i]-In[i])
        return Index

    def Spin2Index(self, In, Out):
        return In*SPIN+Out
    def Spin4Index(self, InTuple, OutTuple):
        return InTuple[IN]*SPIN3+InTuple[OUT]*SPIN2+ \
                OutTuple[IN]*SPIN+OutTuple[OUT]

class Weight:
    def __init__(self, Name, Beta, L, IsSymmetric=True):
        self.__Name=Name
        self.__SmoothTName="{0}.SmoothT".format(Name)
        self.__DeltaTName="{0}.DeltaT".format(Name)
        self.__Beta=Beta
        self.__L=L
        self.__IsSymmetric=IsSymmetric
        self.SmoothT=None
        self.DeltaT=None
    def GetMap(self):
        return IndexMap(self.__Beta, self.__L)

    def fftTime(self,BackForth):
        if BackForth==1:
            self.ChangeSymmetry(1)
            if self.SmoothT is not None:
                self.SmoothT=np.fft.fft(self.SmoothT, axis=TAU)
        if BackForth==-1:
            if self.SmoothT is not None:
                self.SmoothT=np.fft.ifft(self.SmoothT, axis=TAU)
            self.ChangeSymmetry(-1)
    def ChangeSymmetry(self,BackForth):
        ''' the transformation has to be done in continuous tau representation, namely using  
        exp(i*Pi*Tau_n/Beta)(e.g. exp(i*Pi*(n+1/2)/N)) as the phase factor
        otherwise, if you use exp(i*Pi*n/N)) as the phase factor here, you'll have to take care of 
        an extra coeffecient exp(-i*Pi/(2N)) for each function (G0, Sigma, G) in the integral.
        '''
        if self.__IsSymmetric:
            return
        Map=self.GetMap()
        tau=np.array([Map.IndexToTau(e) for e in range(SHAPE[TAU])])
        PhaseFactor=np.exp(1j*BackForth*np.pi*tau/self.__Beta)
        if self.SmoothT is not None:
            self.SmoothT*=PhaseFactor

    def fftSpace(self,BackForth):
        OldShape=self.SmoothT.shape
        Axis, NewShape=self.__SpatialShape(OldShape)
        self.SmoothT.reshape(NewShape)
        if BackForth==1:
            self.SmoothT=np.fft.fftn(self.SmoothT, axis=Axis)   
        elif BackForth==-1:
            self.SmoothT=np.fft.ifftn(self.SmoothT, axis=Axis)   
        self.SmoothT.reshape(OldShape)

        OldShape=self.DeltaT.shape
        Axis, NewShape=self.__SpatialShape(OldShape)
        self.DeltaT.reshape(NewShape)
        if BackForth==1:
            self.DeltaT=np.fft.fftn(self.DeltaT, axis=Axis)   
        elif BackForth==-1:
            self.DeltaT=np.fft.ifftn(self.DeltaT, axis=Axis)   
        self.DeltaT.reshape(OldShape)
    
    def __SpatialShape(self,shape):
        InsertPos=VOL
        SpatialShape=shape[0:InsertPos]+self.__L+shape[InsertPos+1:]
        return range(InsertPos, InsertPos+len(self.__L)), SpatialShape

    def __LoadNpz(self, FileName):
        try:
            data=np.load(FileName)
        except IOError:
            log.error(FileName+" fails to read!")
            sys.exit(0)
        return data

    def Load(self, FileName):
        log.info("Loading {0} Matrix...".format(self.__Name));
        data=self.__LoadNpz(FileName)
        if self.__SmoothTName in data.files:
            log.info("Load {0}".format(self.__SmoothTName))
            self.SmoothT=data[self.__SmoothTName]
            CheckShape(self.SmoothT.shape)
        if self.__DeltaTName in data.files:
            log.info("Load {0}".format(self.__DeltaTName))
            self.DeltaT=data[self.__DeltaTName]
            CheckShape(self.DeltaT.shape)

    def Save(self, FileName, Mode="a"):
        log.info("Saving {0} Matrix...".format(self.__Name));
        data={}
        if Mode is "a" and os.path.exists(FileName)==True:
            olddata=self.__LoadNpz(FileName)
            for e in olddata.files:
                data[e]=olddata[e]
        if self.DeltaT is not None:
            data[self.__DeltaTName]=self.DeltaT
        if self.SmoothT is not None:
            data[self.__SmoothTName]=self.SmoothT
        np.savez(FileName, **data)

class TestWeight(unittest.TestCase):
    def setUp(self):
        self.L=[8,8]
        self.Beta=1.0
        self.G=Weight("G", self.Beta, self.L)
        self.G.SmoothT=np.array([1.0, 2.0])
        self.W=Weight("W", self.Beta, self.L)
        self.W.DeltaT=np.array([2.0, 3.0])

    def test_matrix_IO(self):
        FileName="test.npz"
        self.G.Save(FileName)
        self.W.Save(FileName, "a")
        newW=Weight("W", self.Beta, self.L)
        newW.Load(FileName)
        self.assertTrue(np.allclose(self.W.DeltaT,newW.DeltaT))

class TestWeightFFT(unittest.TestCase):
    def setUp(self):
        self.L=[8,8]
        self.Beta=1.0
        SHAPE[VOL]=self.L[0]*self.L[1]
        SHAPE[TAU]=64
        self.G=Weight("G", self.Beta, self.L, False)
        self.G.SmoothT=np.zeros([4, SHAPE[SUB], SHAPE[VOL], SHAPE[TAU]])+0j
        Map=self.G.GetMap()
        TauGrid=np.linspace(0.0, self.Beta, SHAPE[TAU], endpoint=False)/self.Beta
        #last point<self.Beta!!!
        #print TauGrid
        self.gTau=np.exp(TauGrid)+np.exp(1.0-TauGrid)
        self.G.SmoothT+=self.gTau
    def test_fft_backforth(self):
        self.G.fftTime(1)
        self.G.fftTime(-1)
        self.assertTrue(np.allclose(self.G.SmoothT[0,0,0,:], self.gTau))
    def test_fft_symmetry(self):
        self.G.ChangeSymmetry(-1)
        self.G.fftTime(1) #fftTime(1) will call ChangeSymmetry(1)
        self.assertTrue(np.allclose(self.G.SmoothT[0,0,0,:], np.fft.fft(self.gTau)))
        self.G.fftTime(-1)
        self.G.ChangeSymmetry(1)

if __name__=="__main__":
    G=Weight("G", 1.0,[8,8], False);
    G.Load("../data/GW.npz");
    G.Save("../data/GW_new.npz");
    print G.SmoothT.shape


