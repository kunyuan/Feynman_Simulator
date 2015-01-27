#!/usr/bin/env python
import numpy as np
import sys, os, unittest, math
from logger import *

SPIN,SPIN2,SPIN3=2,4,8
IN,OUT=0,1
DOWN,UP=0,1

class IndexMap:
    def __init__(self, Beta, L, NSublat, MaxTauBin):
        self.MaxTauBin=MaxTauBin
        self.Beta=Beta
        self.Dim=len(L)
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

    def GetAllSublatTuple(self):
        return ((i,j) for i in range(self.NSublat) for j in range(self.NSublat))

    def GetLocalSublatTuple(self):
        return ((i,i) for i in range(self.NSublat))

    def CoordiIndex(self, Out):
        #Out[0]*L1*L2+Out[1]*L2+Out[2] with In=(0,0,0)
        Index=Out[0]
        for i in range(1,len(self.L)):
            Index=Index*self.L[i]+Out[i]
        return Index
    def IndexToCoordi(self, Index):
        Out=[0]*len(self.L)
        for i in range(len(self.L)-1, 0, -1):
            Out[i]=Index%self.L[i]
            Index=Index/self.L[i]
        Out[0]=Index
        return Out

    def GetAllCoordi(self):
        if len(self.L)==2:
            return ((i,j) for i in range(self.L[0]) for j in range(self.L[1]))
        elif len(self.L)==3:
            return ((i,j,k) for i in range(self.L[0]) 
                            for j in range(self.L[1])
                            for k in range(self.L[2]))
        else:
            Assert(False, "Dimension {0} has been implemented!".format(len(self.L)))
    def LatIndex(self, Coordi, Sublat):
        return self.CoordiIndex(Coordi)*self.NSublat+Sublat

    def Spin2Index(self, SpinIn, SpinOut):
        return SpinIn*SPIN+SpinOut

    def Pauli(self):
        return (np.array([[0, 1], [-1, 0]]), np.array([[0, 1j], [-1j, 0]]), np.array([[1,0],[0,-1]]))

    def IsConserved(self, SpinNum, SpinTuple):
        if SpinNum==2:
            return (SpinTuple[IN]==SpinTuple[OUT])
        if SpinNum==4:
            return (SpinTuple[IN][IN]+SpinTuple[OUT][IN]==SpinTuple[IN][OUT]+SpinTuple[OUT][OUT])

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
    """Assumpations of Data:
       TAU axis always follows VOL axis
    """
    def __init__(self, Name, Map, NSpin, Symmetry=None, SpaceDomain="R", TimeDomain="T"):
        """Name: 'SmoothT' or 'DeltaT'
           NSpin: 'NoSpin', 'TwoSpins' or 'FourSpins'
           Symmetry: 'Symmetric' or 'AntiSymmetric', only be checked if TauDep is 'SmoothT'
           SpaceDomain: 'R' or 'K'
           TimeDomain: 'T' or 'W'
        """
        self.Map=Map
        self.Name=Name
        self.NSublat=self.Map.NSublat
        self.L=self.Map.L
        self.Vol=self.Map.Vol
        self.Beta=self.Map.Beta
        self.MaxTauBin=self.Map.MaxTauBin
        if NSpin is "NoSpin":
            self.NSpin=1 #the spin dimension has only one element
        elif NSpin is "TwoSpins":
            self.NSpin=2
        elif NSpin is "FourSpins":
            self.NSpin=4
        else:
            Assert(False, "Only accept NoSpin, TwoSpins or FourSpins, not {0}".format(NSpin))
        self.Shape=[self.NSpin, self.NSublat, self.NSpin, self.NSublat, self.Map.Vol]
        self.VOLDIM=4
        self.TAUDIM=self.VOLDIM+1
        if Name=="SmoothT":
            self.Beta=self.Map.Beta
            self.Shape.append(self.Map.MaxTauBin)
            self.__HasTau=True
            self.__SpaceTimeIndex=[[k,v] for k in range(self.Shape[self.VOLDIM]) 
                                         for v in range(self.Shape[self.TAUDIM])]
            self.__SpaceTimeVol=self.Shape[self.VOLDIM]*self.Shape[self.TAUDIM]

            if Symmetry is "Symmetric":
                self.IsSymmetric=True
            elif Symmetry is "AntiSymmetric":
                self.IsSymmetric=False
            else:
                Assert(False, "Should be either Symmetric or AntiSymmetric, not {0}".format(Symmetry))
        elif Name=="DeltaT":
            self.__HasTau=False
            self.__SpaceTimeIndex=[[k,] for k in range(self.Shape[self.VOLDIM])]
            self.__SpaceTimeVol=self.Shape[self.VOLDIM]
        else:
            Assert(False, "Should be either .SmoothT or .DeltaT in Name, not {0}".format(Name))

        self.Data=np.zeros(self.Shape, dtype=complex)
        self.__OriginShape=list(self.Shape) #get a copy of self.Shape
        Assert(SpaceDomain in ['R', 'K'], "SpaceDomain is either R or K")
        Assert(TimeDomain in ['T', 'W'], "TimeDomain is either T or W")
        self.SpaceDomain=SpaceDomain
        self.TimeDomain=TimeDomain
        
    def Copy(self):
        """return a deep copy of Weight instance"""
        import copy
        return copy.deepcopy(self)

    def FFT(self, *SpaceOrTime):
        if "R" in SpaceOrTime and self.SpaceDomain is "K":
            self.__fftSpace(-1)
            self.SpaceDomain="R"
        if "K" in SpaceOrTime and self.SpaceDomain is "R":
            self.__fftSpace(1)
            self.SpaceDomain="K"
        if "T" in SpaceOrTime and self.TimeDomain is "W":
            self.__fftTime(-1)
            self.TimeDomain="T"
        if "W" in SpaceOrTime and self.TimeDomain is "T":
            self.__fftTime(1)
            self.TimeDomain="W"
    def ChangeSymmetry(self,BackForth):
        ''' the transformation has to be done in continuous tau representation, namely using  
        exp(-i*Pi*Tau_n/Beta)(e.g. exp(-i*Pi*(n+1/2)/N)) as the phase factor
        otherwise, if you use exp(-i*Pi*n/N)) as the phase factor here, you'll have to take care of 
        an extra coeffecient exp(-i*Pi/(2N)) for each function (G0, Sigma, G) in the integral.
        '''
        if self.IsSymmetric or not self.__HasTau:
            return
        tau=np.array([self.Map.IndexToTau(e) for e in range(self.Shape[self.TAUDIM])])
        PhaseFactor=np.exp(-1j*BackForth*np.pi*tau/self.Beta)
        self.Data*=PhaseFactor

    def Inverse(self):
        self.__InverseSpinAndSublat()

    def FromDict(self, data):
        log.info("Loading {0} Matrix...".format(self.Name));
        if self.Name in data:
            log.info("Load {0}".format(self.Name))
            self.Data=data[self.Name]
        else:
            Assert(False, "{0} not found!").format(self.Name)
        return self
    def ToDict(self):
        log.info("Saving {0} Matrix...".format(self.Name));
        return {self.Name: self.Data}

    def __InverseSpinAndSublat(self):
        Sp, Sub = self.NSpin, self.NSublat
        self.Data = self.Data.reshape([Sp*Sub, Sp*Sub, self.__SpaceTimeVol])
        for index in range(self.__SpaceTimeVol):
            try:
                self.Data[:,:,index] = np.linalg.inv(self.Data[:,:,index])
            except:
                log.error("Fail to inverse matrix :,:,{0}\n{1}".format(index, self.Data[index].shape))
                sys.exit(0)
        self.Data = self.Data.reshape(self.__OriginShape)
    def __AssertShape(self, shape1, shape2):
        Assert(tuple(shape1)==tuple(shape2), \
                "Shape {0} is expected instead of shape {1}!".format(shape1, shape2))
    def __fftTime(self,BackForth):
        if not self.__HasTau:
            return
        if BackForth==1:
            self.ChangeSymmetry(1)
            self.Data=np.fft.fft(self.Data, axis=self.TAUDIM)
            self.__AdditionalPhaseFactor(1)
        if BackForth==-1:
            self.__AdditionalPhaseFactor(-1)
            self.Data=np.fft.ifft(self.Data, axis=self.TAUDIM)
            self.ChangeSymmetry(-1)
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
    def __AdditionalPhaseFactor(self,BackForth):
        ''' the transformation has to be done in continuous tau representation, namely using  
        exp(-i*2*Pi*m*Tau_n/Beta)(e.g. exp(-i*2*Pi*m*(n+1/2)/N)) as the phase factor
        so we have to multiply an additional phase factor exp(-i*Pi*m/N)
        '''
        if not self.__HasTau:
            return
        omega=np.array(range(self.Shape[self.TAUDIM]))
        PhaseFactor=np.exp(-1j*BackForth*np.pi*omega/self.Shape[self.TAUDIM])
        self.Data*=PhaseFactor
    def __SpatialShape(self, shape):
        InsertPos=self.VOLDIM
        shape=list(shape)
        SpatialShape=shape[0:InsertPos]+self.L+shape[InsertPos+1:]
        return range(InsertPos, InsertPos+len(self.L)), SpatialShape


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
        self.G=Weight("SmoothT", self.Map, "TwoSpins", "AntiSymmetric","R","T")
        TauGrid=np.linspace(0.0, self.G.Beta, self.G.Shape[self.G.TAUDIM], endpoint=False)/self.G.Beta
        #last point<self.Beta!!!
        self.gTau=np.exp(TauGrid)
        xx,yy=np.meshgrid(range(self.G.L[0]),range(self.G.L[1]))
        zz=np.exp(xx+yy)
        self.z=zz[:,:, np.newaxis]*self.gTau
        self.G.Data+=self.z.reshape(self.G.Shape[self.G.VOLDIM:])
    def test_matrix_IO(self):
        FileName="test.npz"
        self.G.Save(FileName)
        newG=Weight("SmoothT", self.Map, "TwoSpins", "AntiSymmetric", "R","T")
        newG.Load(FileName)
        self.assertTrue(np.allclose(self.G.Data,newG.Data))
        os.system("rm "+FileName)
    def test_fft_backforth(self):
        self.G.FFT("W")
        self.G.FFT("T")
        self.assertTrue(np.allclose(self.G.Data[0,0,:,:], self.z.reshape(self.G.Shape[self.G.VOLDIM:])))
    def test_fft_symmetry(self):
        self.G.ChangeSymmetry(-1)
        self.G.FFT("W") #fftTime(1) will call ChangeSymmetry(1)
        self.assertTrue(np.allclose(self.G.Data[0,0,0,:], np.fft.fft(self.gTau)))
        self.G.FFT("T")
        self.G.ChangeSymmetry(1)
    def test_fft_spatial(self):
        old=self.G.Data.copy()
        self.G.FFT("K")
        zzz=np.fft.fftn(self.z, axes=(0,1))
        self.assertTrue(np.allclose(self.G.Data[0,0,:,:], zzz.reshape(self.G.Shape[self.G.VOLDIM:])))
        self.G.FFT("R")
        self.assertTrue(np.allclose(self.G.Data, old))
