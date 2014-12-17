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
    def __init__(self, Beta, L, NSublat, MaxTauBin):
        self.MaxTauBin=MaxTauBin
        self.Beta=Beta
        self.L=L
        self.Vol=1
        for e in self.L:
            self.Vol*=e
        self.NSublat=NSublat
        self.__dBeta=Beta/self.MaxTauBin
        self.__dBetaInverse=1.0/self.__dBeta
    def GetPara(self):
        return {"L":self.L, "NSublat":self.NSublat, \
                "Beta":self.Beta, "MaxTauBin": self.MaxTauBin}

    def TauIndex(self, In, Out):
        tau=Out-In
        Index=math.floor(tau*self.__dBetaInverse)
        return Index if tau>=0 else Index+self.MaxTauBin

    def IndexToTau(self, Index):
        return (Index+0.5)*self.__dBeta

    def SublatIndex(self, In, Out):
        return In*self.NSublat+Out
    def GetAllSublat(self):
        return ((i,j) for i in range(self.NSublat) for j in range(self.NSublat))
    def GetLocalSublatIndexs(self):
        return (self.SublatIndex(i,i) for i in range(self.NSublat))

    def CoordiIndex(self, Out):
        #Out[0]*L1*L2+Out[1]*L2+Out[2] with In=(0,0,0)
        Index=Out[0]
        for i in range(1,len(Out)):
            Index=Index*self.L[i]+Out[i]
        return Index
    def GetAllCoordi(self):
        if len(self.L)==2:
            return ((i,j) for i in range(self.L[0]) for j in range(self.L[1]))
        elif len(self.L)==3:
            return ((i,j,k) for i in range(self.L[0]) 
                            for j in range(self.L[1])
                            for k in range(self.L[2]))
        else:
            Assert(False, "Dimension {0} has been implemented!".format(len(self.L)))

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

    def IsConserved(self, SpinNum, SpinTuple):
        if SpinNum==2:
            return (SpinTuple[IN]==SpinTuple[OUT])
        if SpinNum==4:
            return (SpinTuple[IN][IN]+SpinTuple[OUT][IN]==SpinTuple[IN][OUT]+SpinTuple[OUT][OUT])

    def GetConservedSpinIndexs(self, SpinNum):
        if SpinNum=="TwoSpins":
            return [self.Spin2Index(*s) for s in self.GetConservedSpinTuple(SpinNum)]
        if SpinNum=="FourSpins":
            return [self.Spin4Index(*s) for s in self.GetConservedSpinTuple(SpinNum)]

    def GetConservedSpinTuple(self, SpinNum):
        if SpinNum=="TwoSpins":
            return [(DOWN,DOWN),(UP,UP)]
        if SpinNum=="FourSpins":
            return self.GetSpin4SimilarTuples((UP,UP),(UP,UP))+ \
                   self.GetSpin4SimilarTuples((UP,UP),(DOWN,DOWN))+ \
                   self.GetSpin4SimilarTuples((DOWN,UP),(UP,DOWN))
    def GetSpin4SimilarTuples(self, InTuple, OutTuple):
        if InTuple==OutTuple:
            if InTuple[IN]==InTuple[OUT]:
                return [((DOWN,DOWN),(DOWN,DOWN)),((UP,UP),(UP,UP))]
            else:
                return [((DOWN,UP),(DOWN,UP)),((UP,DOWN),(UP,DOWN))]
        else:
            if InTuple[IN]==InTuple[OUT] and OutTuple[IN]==OutTuple[OUT]:
                return [((DOWN,DOWN),(UP,UP)),((UP,UP),(DOWN,DOWN))]
            elif InTuple[IN]!=InTuple[OUT] and OutTuple[IN]!=OutTuple[OUT]:
                return [((DOWN,UP),(UP,DOWN)),((UP,DOWN),(DOWN,UP))]
            else:
                return [((DOWN,UP),(UP,UP)),((UP,DOWN),(UP,UP)),
                        ((DOWN,UP),(DOWN,DOWN)),((UP,DOWN),(DOWN,DOWN)),
                        ((UP,UP),(DOWN,UP)),((UP,UP),(UP,DOWN)),
                        ((DOWN,DOWN),(DOWN,UP)),((DOWN,DOWN),(UP,DOWN))]

class Weight():
    def __init__(self, Name, Map, NSpin, Symmetry=None):
        """Name: end with '.SmoothT' or '.DeltaT'
           NSpin: 'OneSpin' or 'TwoSpin'
           Symmetry: 'Symmetric' or 'AntiSymmetric', only be checked if TauDep is 'SmoothT'
        """
        self.Map=Map
        self.Name=Name
        self.NSublat=self.Map.NSublat
        self.L=self.Map.L
        if NSpin is "TwoSpins":
            self.NSpin=2
        elif NSpin is "FourSpins":
            self.NSpin=4
        else:
            Assert(False, "Only accept OneSpin or TwoSpins, not {0}".format(NSpin))

        self.Shape=[self.NSpin**2, self.NSublat**2, self.Map.Vol]
        if ".SmoothT" in Name:
            self.Beta=self.Map.Beta
            self.Shape.append(self.Map.MaxTauBin)
            self.__HasTau=True
            self.__SpaceTimeIndex=[[k,v] for k in range(self.Shape[VOL]) for v in range(self.Shape[TAU])]
            if Symmetry is "Symmetric":
                self.IsSymmetric=True
            elif Symmetry is "AntiSymmetric":
                self.IsSymmetric=False
            else:
                Assert(False, "Should be either Symmetric or AntiSymmetric, not {0}".format(Symmetry))
        elif ".DeltaT" in Name:
            self.__HasTau=False
            self.__SpaceTimeIndex=[[k,] for k in range(self.Shape[VOL])]
        else:
            Assert(False, "Should be either .SmoothT or .DeltaT in Name, not {0}".format(Name))

        self.Data=np.zeros(self.Shape, dtype=complex)
        self.__OriginShape=list(self.Shape) #get a copy of self.Shape

    def FFT(self, BackForth, *SpaceOrTime):
        if "Space" in SpaceOrTime:
            self.__fftSpace(BackForth)
        if "Time" in SpaceOrTime:
            self.__fftTime(BackForth)

    def __fftTime(self,BackForth):
        if not self.__HasTau:
            return
        if BackForth==1:
            self.ChangeSymmetry(1)
            self.Data=np.fft.fft(self.Data, axis=TAU)
        if BackForth==-1:
            self.Data=np.fft.ifft(self.Data, axis=TAU)
            self.ChangeSymmetry(-1)
    def ChangeSymmetry(self,BackForth):
        ''' the transformation has to be done in continuous tau representation, namely using  
        exp(-i*Pi*Tau_n/Beta)(e.g. exp(-i*Pi*(n+1/2)/N)) as the phase factor
        otherwise, if you use exp(-i*Pi*n/N)) as the phase factor here, you'll have to take care of 
        an extra coeffecient exp(-i*Pi/(2N)) for each function (G0, Sigma, G) in the integral.
        '''
        if self.IsSymmetric or not self.__HasTau:
            return
        tau=np.array([self.Map.IndexToTau(e) for e in range(self.Shape[TAU])])
        PhaseFactor=np.exp(-1j*BackForth*np.pi*tau/self.Beta)
        self.Data*=PhaseFactor

    def __fftSpace(self, BackForth):
        OldShape=self.__OriginShape
        self.__AssertShape(self.Shape, OldShape)
        Axis, NewShape=self.__SpatialShape(OldShape)
        self.Data=self.Data.reshape(NewShape)
        if BackForth==1:
            self.Data=np.fft.fftn(self.Data, axes=Axis)   
        elif BackForth==-1:
            self.Data=np.fft.ifftn(self.Data, axes=Axis)   
        self.Data=self.Data.reshape(OldShape)
    def __SpatialShape(self, shape):
        InsertPos=VOL
        shape=list(shape)
        SpatialShape=shape[0:InsertPos]+self.L+shape[InsertPos+1:]
        return range(InsertPos, InsertPos+len(self.L)), SpatialShape

    def Reshape(self, ShapeMode):
        OriginShape=self.__OriginShape
        MidShape1=[self.NSpin, self.NSpin,self.NSublat,self.NSublat]+self.__OriginShape[VOL:]
        MidShape2=[self.NSpin, self.NSublat, self.NSpin, self.NSublat]+self.__OriginShape[VOL:]
        NewShape=[self.NSpin*self.NSublat, self.NSpin*self.NSublat]+self.__OriginShape[VOL:]
        if ShapeMode is "SP2SUB2":
            self.__AssertShape(self.Shape, NewShape)
            self.Data=self.Data.reshape(MidShape2).swapaxes(1,2).reshape(OriginShape)
            self.Shape=OriginShape
        elif ShapeMode is "SPSUBSPSUB":
            self.__AssertShape(self.Shape, OriginShape)
            self.Data=self.Data.reshape(MidShape1).swapaxes(1,2).reshape(NewShape)
            self.Shape=NewShape

    def Inverse(self):
        self.__InverseSpinAndSublat()

    def __InverseSpinAndSublat(self):
        for j in self.__SpaceTimeIndex:
            index=[Ellipsis,]+j
            try:
                self.Data[index] = np.linalg.inv(self.Data[index])
            except:
                log.error("Fail to inverse matrix :,:,{0}\n{1}".format(index, self.Data[index].shape))
                sys.exit(0)

    def Load(self, FileName):
        log.info("Loading {0} Matrix...".format(self.Name));
        data=self.__LoadNpz(FileName)
        if self.Name in data.files:
            log.info("Load {0}".format(self.Name))
            self.__AssertShape(self.Shape, data[self.Name].shape)
            self.Data=data[self.Name]
        else:
            Assert(False, "{0} not found!").format(self.Name)
    def Save(self, FileName, Mode="a"):
        log.info("Saving {0} Matrix...".format(self.Name));
        data={}
        if Mode is "a" and os.path.exists(FileName)==True:
            olddata=self.__LoadNpz(FileName)
            for e in olddata.files:
                data[e]=olddata[e]
        data[self.Name]=self.Data
        np.savez(FileName, **data)
    def __LoadNpz(self, FileName):
        try:
            data=np.load(FileName)
        except IOError:
            log.error(FileName+" fails to read!")
            sys.exit(0)
        return data
    def __AssertShape(self, shape1, shape2):
        Assert(tuple(shape1)==tuple(shape2), \
                "Shape {0} is expected instead of shape {1}!".format(shape1, shape2))

class TestIndexMap(unittest.TestCase):
    def setUp(self):
        self.L=[8,8]
        self.Beta=1.0
        self.Map=IndexMap(self.Beta, self.L, NSublat=2, MaxTauBin=64)
    def test_conserved_spin_filter(self):
        for s in self.Map.GetConservedSpinTuple("TwoSpins"):
            self.assertTrue(s[IN]==s[OUT])
        for s in self.Map.GetConservedSpinTuple("FourSpins"):
            self.assertTrue(s[IN][IN]+s[OUT][IN]==s[OUT][OUT]+s[IN][OUT])

class TestWeightFFT(unittest.TestCase):
    def setUp(self):
        self.Map=IndexMap(Beta=1.0, L=[8,8], NSublat=2, MaxTauBin=64)
        self.G=Weight("G.SmoothT", self.Map, "TwoSpins", "AntiSymmetric")
        TauGrid=np.linspace(0.0, self.G.Beta, self.G.Shape[TAU], endpoint=False)/self.G.Beta
        #last point<self.Beta!!!
        self.gTau=np.exp(TauGrid)
        xx,yy=np.meshgrid(range(self.G.L[0]),range(self.G.L[1]))
        zz=np.exp(xx+yy)
        self.z=zz[:,:, np.newaxis]*self.gTau
        self.G.Data+=self.z.reshape(self.G.Shape[VOL:])
    def test_matrix_IO(self):
        FileName="test.npz"
        self.G.Save(FileName)
        newG=Weight("G.SmoothT", self.Map, "TwoSpins", "AntiSymmetric")
        newG.Load(FileName)
        self.assertTrue(np.allclose(self.G.Data,newG.Data))
        os.system("rm "+FileName)
    def test_fft_backforth(self):
        self.G.FFT(1, "Time")
        self.G.FFT(-1, "Time")
        self.assertTrue(np.allclose(self.G.Data[0,0,:,:], self.z.reshape(self.G.Shape[VOL:])))
    def test_fft_symmetry(self):
        self.G.ChangeSymmetry(-1)
        self.G.FFT(1, "Time") #fftTime(1) will call ChangeSymmetry(1)
        self.assertTrue(np.allclose(self.G.Data[0,0,0,:], np.fft.fft(self.gTau)))
        self.G.FFT(-1, "Time")
        self.G.ChangeSymmetry(1)
    def test_fft_spatial(self):
        old=self.G.Data.copy()
        self.G.FFT(1, "Space")
        zzz=np.fft.fftn(self.z, axes=(0,1))
        self.assertTrue(np.allclose(self.G.Data[0,0,:,:], zzz.reshape(self.G.Shape[VOL:])))
        self.G.FFT(-1, "Space")
        self.assertTrue(np.allclose(self.G.Data, old))
