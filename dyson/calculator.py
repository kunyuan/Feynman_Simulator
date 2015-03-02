#!usr/bin/env python
import sys
import numpy as np
import parameter as para
from weight import UP,DOWN,IN,OUT
import weight, plot
from logger import *
print "calculator"

def SigmaSmoothT_FirstOrder(G, W, map):
    '''Fock diagram, assume Spin Conservation'''
    OrderSign=-1
    Sigma=weight.Weight("SmoothT", map, "TwoSpins", "AntiSymmetric", "R", "T")
    TauRange = range(map.MaxTauBin)
    G.FFT("R", "T")
    W.FFT("R", "T")
    for spin1 in range(2):
        for spin2 in range(2):
            spinW = (map.Spin2Index(spin1,spin2),map.Spin2Index(spin2,spin1))
            spinG = (spin2, spin2)
            spinSigma = (spin1, spin1)
            Sigma.Data[spinSigma[IN], :, spinSigma[OUT], :, :]  \
                    += OrderSign*G.Data[spinG[IN], :, spinG[OUT], :, :]\
                    *W.Data[spinW[IN], :, spinW[OUT], :, :]
    return Sigma

def SigmaDeltaT_FirstOrder(G, W0, map):
    '''Hatree-Fock diagram, assume Spin Conservation'''
    ########Fock Diagram
    OrderSign=-1
    AntiSymmetricFactor=-1
    SigmaDeltaT=weight.Weight("DeltaT", map, "TwoSpins", "AntiSymmetric", "R")
    G.FFT("R", "T")
    W0.FFT("R")
    for spin1 in range(2):
        for spin2 in range(2):
            #############G(tau==-0)
            spinW = (map.Spin2Index(spin1,spin2),map.Spin2Index(spin2,spin1))
            spinG = (spin2, spin2)
            spinSigma = (spin1, spin1)
            SigmaDeltaT.Data[spinSigma[IN], :, spinSigma[OUT], :, :]+= OrderSign \
                    *AntiSymmetricFactor*G.Data[spinG[IN], :, spinG[OUT], :, :, -1]\
                    *W0.Data[spinW[IN], :, spinW[OUT], :, :]

    ########Hatree Diagram, or bubble diagram
    FermiLoopSign=-1
    for sp1 in range(2):
        for sp2 in range(2):
            spinW = (map.Spin2Index(sp1,sp1), map.Spin2Index(sp2,sp2))
            for sub1 in range(map.NSublat):
                for sub2 in range(map.NSublat):
                    for r in range(map.Vol):
                        SigmaDeltaT.Data[sp1, sub1, sp1, sub1, 0]+= OrderSign*FermiLoopSign \
                            *AntiSymmetricFactor*G.Data[sp2, sub2, sp2, sub2, 0, -1] \
                            *W0.Data[spinW[IN], sub1, spinW[OUT], sub2, r]
    return SigmaDeltaT


def Polar_FirstOrder(G, map):
    OrderSign=-1
    FermiLoopSign=-1
    AntiSymmetricFactor=-1
    Polar=weight.Weight("SmoothT", map, "FourSpins","Symmetric", "R","T")
    G.FFT("R","T")
    NSublat = map.NSublat
    NSpin = G.NSpin
    SubList=[(a,b) for a in range(NSublat) for b in range(NSublat)]
    SpList=[(a,b) for a in range(NSpin) for b in range(NSpin)]
    for spin1,spin2 in SpList:
        spinPolar = ((spin1,spin2),(spin2,spin1))
        spinG1 = (spin1, spin1)
        spinG2 = (spin2, spin2)
        for subA,subB in SubList:
            Polar.Data[map.Spin2Index(*spinPolar[IN]),subA, \
                    map.Spin2Index(*spinPolar[OUT]),subB,:,:]+=OrderSign*FermiLoopSign \
                    *AntiSymmetricFactor*G.Data[spinG1[IN], subB, spinG1[OUT], subA, :, ::-1]  \
                    *G.Data[spinG2[IN], subA, spinG2[OUT], subB, :, :]
    return Polar

def W_FirstOrder(W0, Polar, map):
    W=weight.Weight("SmoothT", map, "FourSpins", "Symmetric", "K","W")
    W0.FFT("K")
    Polar.FFT("K","W")
    Beta=map.Beta
    TauRange = range(map.MaxTauBin)
    SubRange=range(W.NSublat)
    SubList=[(a,b,c,d) for a in SubRange for b in SubRange for c in SubRange for d in SubRange]
    SpinList=[(Wtuple,Polartuple) for Wtuple in map.GetConservedSpinTuple("FourSpins") \
                                  for Polartuple in map.GetConservedSpinTuple("FourSpins") \
                                  if map.IsConserved(4, (Wtuple[IN], Polartuple[IN]))] 
    #make sure spin conservation on W0

    for spWt,spPolart in SpinList:
        spW0L=(map.Spin2Index(*spWt[IN]), map.Spin2Index(*spPolart[IN]))
        spW0R=(map.Spin2Index(*spPolart[OUT]), map.Spin2Index(*spWt[IN]))
        spW = (map.Spin2Index(*spWt[IN]), map.Spin2Index(*spWt[OUT]))
        spPolar = (map.Spin2Index(*spPolart[IN]), map.Spin2Index(*spPolart[OUT]))
        for e in SubList:
            for tau in TauRange:
                W.Data[spW[IN],e[0],spW[OUT],e[3],:,tau]+=\
                        W0.Data[spW0L[IN],e[0],spW0L[OUT],e[1],:] \
                        *Polar.Data[spPolar[IN],e[1],spPolar[OUT],e[2],:,tau]\
                        *W0.Data[spW0R[IN],e[2],spW0R[OUT],e[3],:]
    return W

def Calculate_Denorminator(W0, Polar, map):
    """require W0 and Polar to be in k, omega space"""
    W0.FFT("K")
    Polar.FFT("K","W")
    Beta=map.Beta
    NSpin, NSub=W0.NSpin, W0.NSublat
    JP=np.einsum("ijklv,klmnvt->ijmnvt",W0.Data, Polar.Data)
    #JP shape: NSpin,NSub,NSpin,NSub,Vol,Tau
    Temp=JP*Beta/map.MaxTauBin
    for tau in range(map.MaxTauBin):
        Temp[...,tau] *= np.cos(tau*np.pi/map.MaxTauBin)
    I=np.eye(NSpin*NSub).reshape([NSpin,NSub,NSpin,NSub])
    return I[...,np.newaxis,np.newaxis]-Temp, JP

def W_Dyson(W0, Polar, map, Lat):
    W=weight.Weight("SmoothT", map, "FourSpins", "Symmetric", "K","W")
    ChiTensor=weight.Weight("SmoothT", map, "FourSpins", "Symmetric", "K","W")

    Polar.FFT("K","W")
    Denorm,JP=Calculate_Denorminator(W0, Polar, map)

    W0.FFT("K")
    JPJ=np.einsum("ijklvt,klmnv->ijmnvt", JP, W0.Data)
    lu_piv,Determ=weight.LUFactor(Denorm)
    Check_Denorminator(Determ,map)
    W.LUSolve(lu_piv, JPJ)
    ChiTensor.LUSolve(lu_piv, -Polar.Data)
    return W, ChiTensor, Determ

def G_Dyson(G0, SigmaDeltaT, Sigma, map):
    Beta=map.Beta
    G=weight.Weight("SmoothT", map, "TwoSpins", "AntiSymmetric", "K","W")
    G0.FFT("K", "W")
    SigmaDeltaT.FFT("K")
    Sigma.FFT("K", "W")

    NSpin, NSub=G.NSpin, G.NSublat

    G0SigmaDeltaT=np.einsum("ijklvt,klmnv->ijmnvt",G0.Data, SigmaDeltaT.Data)
    G0Sigma=np.einsum("ijklvt,klmnvt->ijmnvt",G0.Data, Sigma.Data)

    ####correction term
    for tau in range(map.MaxTauBin):
        G0SigmaDeltaT[...,tau]*= np.cos(np.pi*map.IndexToTau(tau)/Beta)

    GS  = Beta/map.MaxTauBin*(Beta/map.MaxTauBin*G0Sigma+G0SigmaDeltaT) 
    #GS shape: NSpin,NSub,NSpin,NSub,Vol,Tau

    I=np.eye(NSpin*NSub).reshape([NSpin,NSub,NSpin,NSub])
    Denorm=I[...,np.newaxis,np.newaxis]-GS
    lu_piv,Determ=weight.LUFactor(Denorm)
    Check_Denorminator(Determ,map)
    G.LUSolve(lu_piv, G0.Data);
    return G

def Add_ChiTensor_ZerothOrder(ChiTensor, G, map):
    """add G(tau=0^-, k=0)G(tau=0^-, k=0) to ChiTensor
       This diagram has +1 sign since two fermi loop contribute (-1)^2
    """
    G.FFT("R", "T")
    ChiTensor.FFT("R", "T")
    NSub=map.NSublat
    SpList=[(a,b,c,d) for a in range(2) for b in range(2) for c in range(2) for d in range(2)]
    SubList=[(a,b) for a in range(NSub) for b in range(NSub)]
    ChiTensorInUnitCell=np.empty((ChiTensor.NSpin,NSub,ChiTensor.NSpin,NSub), dtype=complex)
    for spIn1, spIn2, spOut1, spOut2 in SpList:
        for subIn, subOut in SubList:
            for tau in range(map.MaxTauBin):
                ChiSpIn=map.Spin2Index(spIn1,spIn2)
                ChiSpOut=map.Spin2Index(spOut1,spOut2)
                ChiTensorInUnitCell[ChiSpIn,subIn,ChiSpOut,subOut]=G.Data[spIn1,subIn,spIn2,subIn,0,-1]*G.Data[spOut1,subOut,spOut2,subOut,0,-1]
    ChiTensor.Data+=ChiTensorInUnitCell[...,np.newaxis,np.newaxis]
    return ChiTensor

def Calculate_Chi(ChiTensor, map):
    NSpin, NSublat=ChiTensor.NSpin, ChiTensor.NSublat
    SxSx=np.zeros((NSpin,NSpin))
    SySy=np.zeros((NSpin,NSpin))
    SzSz=np.zeros((NSpin,NSpin))
    UU=map.Spin2Index(UP,UP) 
    UD=map.Spin2Index(UP,DOWN) 
    DU=map.Spin2Index(DOWN, UP) 
    DD=map.Spin2Index(DOWN, DOWN) 
    SxSx[UD, UD]=SxSx[DU, DU]=1
    SxSx[UD, DU]=SxSx[DU, UD]=1
    SySy[UD, UD]= SySy[DU, DU]=-1
    SySy[UD, DU]= SySy[DU, UD]=1
    SzSz[UU, UU]= SzSz[DD, DD]=1
    SzSz[UU, DD]= SzSz[DD, UU]=-1
    Chi=weight.Weight("SmoothT", map, "NoSpin", "Symmetric", ChiTensor.SpaceDomain, ChiTensor.TimeDomain)
    Chi.Data=0.0
    SS=[SxSx/4.0, SySy/4.0, SzSz/4.0]
    for i in range(3):
        temp=np.einsum("ik, kminvt->mnvt", SS[i], ChiTensor.Data)
        Chi.Data+=temp.reshape([1, NSublat, 1, NSublat, map.Vol, map.MaxTauBin]) 
    return Chi

class DenorminatorTouchZero(Exception):
    def __init__(self, value, pos, freq):
       self.value = value
       self.position=pos
       self.frequency=freq
    def __str__(self):
        return "Denorminator touch zero, minmum {0} is at K={1} and Omega={2}".format(self.value, self.position, self.frequency)

def Check_Denorminator(Determ, map):
    pos=np.where(Determ==Determ.min())
    x,t=pos[0][0], pos[1][0]
    log.info("The minmum {0} is at K={1} and Omega={2}".format(Determ.min(), map.IndexToCoordi(x), t))
    if Determ.min().real<0.0 and Determ.min().imag<1.0e-4:
        raise DenorminatorTouchZero(Determ.min(), map.IndexToCoordi(x), t)
