#!/usr/bin/env python
import math
import numpy as np
import sys
import os
import unittest
from logger import *

SPIN,SPIN2,SPIN3=2,4,8
IN,OUT=0,1
DOWN,UP=0,1
SP,SUB,VOL,TAU=0,1,2,3

### Shape: [SP][SUB][VOL][TAU], the TAU dimension may be missing

class IndexMap:
    def __init__(self, Beta, L_, NSublattice, MAX_TAU_BIN):
        self.__MaxTauBin=MAX_TAU_BIN
        self.__Beta=Beta
        self.__dBeta=Beta/self.__MaxTauBin
        self.__dBetaInverse=1.0/self.__dBeta
        self.__L=L_
        self.__NSublattice=NSublattice

    def TauIndex(self, In, Out):
        tau=Out-In
        Index=math.floor(tau*self.__dBetaInverse)
        return Index if tau>=0 else Index+self.__MaxTauBin

    def IndexToTau(self, Index):
        return (Index+0.5)*self.__dBeta

    def SublatIndex(self, In, Out):
        return In*self.__NSublattice+Out

    def CoordiIndex(self,In, Out):
        #Out[0]*L1*L2+Out[1]*L2+Out[2] with In=(0,0,0)
        Index=Out[0]-In[0]
        for i in range(1,len(Vec)):
            Index=Index*self.__L[i]+(Out[i]-In[i])
        return Index

    def SpinIndex(self, In, Out, SpinNum):
        if SpinNum==2:
            self.Spin2Index(In, Out)
        elif SpinNum==4:
            self.Spin4Index(In, Out)
    def Spin2Index(self, In, Out):
        return In*SPIN+Out
    def Spin4Index(self, InTuple, OutTuple):
        return InTuple[IN]*SPIN3+InTuple[OUT]*SPIN2+ \
                OutTuple[IN]*SPIN+OutTuple[OUT]

    def IsLegal(self, SpinNum, SpinTuple):
        if SpinNum==2:
            return (SpinTuple[IN]==SpinTuple[OUT])
        if SpinNum==4:
            return (SpinTuple[IN][IN]+SpinTuple[OUT][IN]==SpinTuple[IN][OUT]+SpinTuple[OUT][OUT])

    def GetLegalSpinIndexs(self, SpinNum):
        if SpinNum==2:
            return [self.Spin2Index(*s) for s in self.GetLegalSpinTuple(SpinNum)]
        if SpinNum==4:
            return [self.Spin4Index(*s) for s in self.GetLegalSpinTuple(SpinNum)]

    def GetLegalSpinTuple(self, SpinNum):
        if SpinNum==2:
            return [(DOWN,DOWN),(UP,UP)]
        if SpinNum==4:
            return [((DOWN,DOWN),(DOWN,DOWN)),((UP,UP),(UP,UP)), \
                    ((UP,UP),(DOWN,DOWN)),((DOWN,DOWN),(UP,UP)), \
                    ((UP,DOWN),(DOWN,UP)),((DOWN,UP),(UP,DOWN))]

class Weight:
    def __init__(self, Name, Beta, L, IsSymmetric=True):
        self.__Name=Name
        self.__SmoothTName="{0}.SmoothT".format(Name)
        self.__DeltaTName="{0}.DeltaT".format(Name)
        self.__Beta=Beta
        self.__L=L
        self.__IsSymmetric=IsSymmetric
        self.__Shape=[None,]*4
        self.SmoothT=None
        self.DeltaT=None
    def SetShape(self, shape_):
        for e in range(len(self.__Shape)):
            if e>=len(shape_):
                pass
            elif self.__Shape[e] is None:
                self.__Shape[e]=shape_[e]
            else:
                Assert(self.__Shape[e]==shape_[e], \
                        "Shape {0} does not match {1} at the {2}rd element!".format(self.__Shape, shape_,e))
        self.__SpinNum=int(math.sqrt(self.__Shape[SP]))
        self.__NSublattice=int(math.sqrt(self.__Shape[SUB]))
    def GetMap(self):
        return IndexMap(self.__Beta, self.__L, self.__NSublattice, self.__Shape[TAU])

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
        tau=np.array([Map.IndexToTau(e) for e in range(self.__Shape[TAU])])
        PhaseFactor=np.exp(1j*BackForth*np.pi*tau/self.__Beta)
        if self.SmoothT is not None:
            self.SmoothT*=PhaseFactor

    def fftSpace(self,BackForth):
        self.SmoothT=self.__fftSpace(self.SmoothT, BackForth)
        self.DeltaT=self.__fftSpace(self.DeltaT, BackForth)
    def __fftSpace(self, array, BackForth):
        if array is None:
            return None
        OldShape=array.shape
        Axis, NewShape=self.__SpatialShape(OldShape)
        temp=array.reshape(NewShape)
        if BackForth==1:
            temp=np.fft.fftn(temp, axes=Axis)   
        elif BackForth==-1:
            temp=np.fft.ifftn(temp, axes=Axis)   
        return temp.reshape(OldShape)
    def __SpatialShape(self, shape):
        InsertPos=VOL
        shape=list(shape)
        SpatialShape=shape[0:InsertPos]+self.__L+shape[InsertPos+1:]
        return range(InsertPos, InsertPos+len(self.__L)), SpatialShape

    def Inverse(self):
        if self.__SpinNum==2:
            self.SmoothT=self.__InverseSublat(self.SmoothT)
            self.DeltaT=self.__InverseSublat(self.DeltaT)
        elif self.__SpinNum==4:
            self.SmoothT=self.__InverseSpinAndSublat(self.SmoothT)
            self.DeltaT=self.__InverseSpinAndSublat(self.DeltaT)

    def __InverseSublat(self, array):
        OldShape=array.shape
        NSublat=self.__NSublattice
        temp=array.reshape(OldShape[SP],NSublat,NSublat,OldShape[VOL]*OldShape[TAU])
        for i in self.GetMap().GetLegalSpinIndexs(self.__SpinNum):
            for j in range(temp.shape[-1]):
                try:
                    temp[i,:,:,j] = np.linalg.inv(temp[i,:,:,j])
                except:
                    log.error("Fail to inverse matrix {0},:,:,{1}\n{2}".format(i,j, temp[i,:,:,j]))
        return temp.reshape(OldShape)

    def __InverseSpinAndSublat(self, array):
        OldShape=array.shape
        NSublat=self.__NSublattice
        SpinNum=self.__SpinNum
        temp=array.reshape(SpinNum,SpinNum,NSublat,NSublat,OldShape[VOL]*OldShape[TAU])
        for j in range(temp.shape[-1]):
            try:
                temp[:,:,:,:,j] = np.linalg.tensorinv(temp[:,:,:,:,j])
            except:
                log.error("Fail to inverse matrix :,:,:,:,{1}\n{2}".format(j, temp[:,:,:,:,j]))
        return temp.reshape(OldShape)

    def Load(self, FileName):
        log.info("Loading {0} Matrix...".format(self.__Name));
        data=self.__LoadNpz(FileName)
        if self.__SmoothTName in data.files:
            log.info("Load {0}".format(self.__SmoothTName))
            self.SmoothT=data[self.__SmoothTName]
            self.SetShape(self.SmoothT.shape)
        if self.__DeltaTName in data.files:
            log.info("Load {0}".format(self.__DeltaTName))
            self.DeltaT=data[self.__DeltaTName]
            self.SetShape(self.DeltaT.shape)
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
    def __LoadNpz(self, FileName):
        try:
            data=np.load(FileName)
        except IOError:
            log.error(FileName+" fails to read!")
            sys.exit(0)
        return data

class TestIndexMap(unittest.TestCase):
    def setUp(self):
        self.L=[8,8]
        self.Beta=1.0
        self.shape=[4, 4, self.L[0]*self.L[1], 64]
        self.Map=IndexMap(self.Beta, self.L, self.shape)

    def test_legal_spin_filter(self):
        for s in self.Map.GetLegalSpinTuple(2):
            self.assertTrue(s[IN]==s[OUT])
        for s in self.Map.GetLegalSpinTuple(4):
            self.assertTrue(s[IN][IN]+s[OUT][IN]==s[OUT][OUT]+s[IN][OUT])

class TestWeightFFT(unittest.TestCase):
    def setUp(self):
        self.L=[8,8]
        self.Beta=1.0
        self.G=Weight("G", self.Beta, self.L, False)
        self.shape=[4, 4, self.L[0]*self.L[1], 64]
        self.G.SmoothT=np.zeros(self.shape)+0j
        self.G.SetShape(self.shape) 
        Map=self.G.GetMap()
        TauGrid=np.linspace(0.0, self.Beta, self.shape[TAU], endpoint=False)/self.Beta
        #last point<self.Beta!!!
        self.gTau=np.exp(TauGrid)
        xx,yy=np.meshgrid(range(self.L[0]),range(self.L[1]))
        zz=np.exp(xx+yy)
        self.z=zz[:,:, np.newaxis]*self.gTau
        self.G.SmoothT+=self.z.reshape(self.shape[VOL:])
    def test_matrix_IO(self):
        FileName="test.npz"
        self.G.Save(FileName)
        newG=Weight("G", self.Beta, self.L, False)
        newG.Load(FileName)
        self.assertTrue(np.allclose(self.G.SmoothT,newG.SmoothT))
        os.system("rm "+FileName)
    def test_fft_backforth(self):
        self.G.fftTime(1)
        self.G.fftTime(-1)
        self.assertTrue(np.allclose(self.G.SmoothT[0,0,:,:], self.z.reshape(self.shape[VOL:])))
    def test_fft_symmetry(self):
        self.G.ChangeSymmetry(-1)
        self.G.fftTime(1) #fftTime(1) will call ChangeSymmetry(1)
        self.assertTrue(np.allclose(self.G.SmoothT[0,0,0,:], np.fft.fft(self.gTau)))
        self.G.fftTime(-1)
        self.G.ChangeSymmetry(1)
    def test_fft_spatial(self):
        self.G.fftSpace(1)
        zzz=np.fft.fftn(self.z, axes=(0,1))
        self.assertTrue(np.allclose(self.G.SmoothT[0,0,:,:], zzz.reshape(self.shape[VOL:])))
        self.G.fftSpace(-1)

if __name__=="__main__":
    G=Weight("G", 1.0,[8,8], False);
    G.Load("../data/GW.npz");
    G.fftSpace(1)
    G.Inverse()
    G.Save("../data/GW_new.npz");
    print G.SmoothT.shape


