#!/usr/bin/env python
import math
import numpy as np
import sys
import os.path
import unittest
from logger import *
MAX_TAU_BIN=128
SublatVol=2
D=2
SPIN=2
SPIN2=4
SPIN3=8
IN=0
OUT=1
DOWN=0
UP=1

class IndexMap:
    def __init__(self,Beta, L):
        self.__Beta=Beta
        self.__dBeta=1.0/Beta
        self.__dBetaInverse=1.0/self.__dBeta
        self.__L=L

    def TauIndex(self, In, Out):
        Index=math.floor((Out-In)*self.__dBetaInverse)
        return Index if tau>=0 else Index+MAX_TAU_BIN
    def IndexToTau(self, Index):
        return (Index+0.5)*self.__dBeta

    def SublatIndex(self, In, Out):
        return In*SublatVol+Out
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
    def __init__(self, Name, Beta, IsSymmetric=True):
        self.__Name=Name
        self.__SmoothTName="{0}.SmoothT".format(Name)
        self.__DeltaTName="{0}.DeltaT".format(Name)
        self.__Beta=Beta
        self.__IsSymmetric=IsSymmetric
        self.SmoothT=None
        self.DeltaT=None
    
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
        if self.__DeltaTName in data.files:
            log.info("Load {0}".format(self.__DeltaTName))
            self.DeltaT=data[self.__DeltaTName]
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
        self.G=Weight("G", 1.0)
        self.G.SmoothT=np.array([1.0, 2.0])
        self.W=Weight("W", 1.0)
        self.W.DeltaT=np.array([2.0, 3.0])
    def test_matrix_IO(self):
        FileName="test.npz"
        self.G.Save(FileName)
        self.W.Save(FileName, "a")
        newW=Weight("W",1.0)
        newW.Load(FileName)
        self.assertTrue(np.allclose(self.W.DeltaT,newW.DeltaT))

if __name__=="__main__":
    G=Weight("G", 1.0, False);
    G.Load("../data/GW.npz");
    G.Save("../data/GW_new.npz");
    print G.SmoothT.shape


