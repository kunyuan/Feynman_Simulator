#!usr/bin/env python
import sys
import numpy as np
import parameter as para
from weight import UP,DOWN,IN,OUT,TAU,SP,SUB,VOL
import weight
from logger import *

def PlotTime(array, Beta):
    import matplotlib.pyplot as plt
    x=np.linspace(0, Beta, len(array))
    plt.plot(x,array,'-')
    plt.show()


def Polar_FirstOrder(G, map):
    Polar=weight.Weight("Polar.SmoothT", map, "FourSpins","Symmetric")
    NSublat = G.NSublat
    SubList=[(map.SublatIndex(a,b),map.SublatIndex(b,a)) for a in range(NSublat) for b in range(NSublat)]
    for spin1 in range(2):
        for spin2 in range(2):
            spinPolar = map.Spin4Index((spin1,spin2),(spin2,spin1))
            spinG1 = map.Spin2Index(spin1, spin1)
            spinG2 = map.Spin2Index(spin2, spin2)
            for subA2B,subB2A in SubList:
                Polar.Data[spinPolar, subA2B, :, :]  \
                        = (-1.0)*G.Data[spinG1, subB2A, :, ::-1]\
                        *G.Data[spinG2, subA2B, :, :]
    return Polar

def W_Dyson(W0,Polar,map):
    W=weight.Weight("W.SmoothT", map, "FourSpins", "Symmetric")
    W0.FFT(1, "Space")
    Polar.FFT(1, "Space", "Time")

    NSpin, NSub=W.NSpin, W.NSublat
    W0.Reshape("SPSUB,SPSUB")
    W.Reshape("SPSUB,SPSUB")
    Polar.Reshape("SPSUB,SPSUB")
    JP=np.einsum("ijv,jkvt->ikvt",W0.Data, Polar.Data)
    #JP shape: NSpin*NSub,NSpin*NSub,Vol,Tau
    I=np.eye(NSpin*NSub)
    W.Data=I[...,np.newaxis,np.newaxis]-JP
    W.Inverse();
    W.Data=np.einsum('ijvt,jkv->ikvt', W.Data,W0.Data)-W0.Data[...,np.newaxis]
    W.Reshape("SP2,SUB2")
    W0.Reshape("SP2,SUB2")
    Polar.Reshape("SP2,SUB2")
    W.FFT(-1, "Space", "Time")
    Polar.FFT(-1, "Space", "Time")
    W0.FFT(-1, "Space")
    #print W.Data[map.Spin4Index((0,0),(0,0)),map.SublatIndex(0,0),0,:]

def Sigma_FirstOrder(G0, W, map):
    Sigma=weight.Weight("Sigma.SmoothT", map, "TwoSpins", "AntiSymmetric")
    return Sigma

def W_FirstOrder(Beta,W0, Polar, map):
    W=weight.Weight("W.SmoothT", map, "FourSpins", "Symmetric")
    TauRange = range(W.Shape[TAU])
    SubRange=range(W.NSublat)
    SubList=[(a,b,c,d) for a in SubRange for b in SubRange for c in SubRange for d in SubRange]
    SpinList=[(Wtuple,Polartuple) for Wtuple in map.GetConservedSpinTuple("FourSpins") \
                                  for Polartuple in map.GetConservedSpinTuple("FourSpins") \
                                  if map.IsConserved(4, (Wtuple[IN], Polartuple[IN]))] 
    #make sure spin conservation on W0

    W0.FFT(1, "Space")
    Polar.FFT(1, "Space")
    for spWt,spPolart in SpinList:
        spW0L=map.Spin4Index(spWt[IN], spPolart[IN])
        spW0R=map.Spin4Index(spPolart[OUT], spWt[IN])
        spW = map.Spin4Index(*spWt)
        spPolar = map.Spin4Index(*spPolart)
        for e in SubList:
            subW0L = map.SublatIndex(e[0], e[1])
            subPolar = map.SublatIndex(e[1], e[2])
            subW0R = map.SublatIndex(e[2], e[3])
            subW = map.SublatIndex(e[0], e[3])
            for tau in TauRange:
                W.Data[spW,subW,:,tau]+=W0.Data[spW0L,subW0L,:] \
                    *Polar.Data[spPolar,subPolar,:,tau]*W0.Data[spW0R,subW0R,:]
    
    W.Data[:,:,:,:] *= (Beta/(W.Shape[TAU]**2.0))
    W0.FFT(-1, "Space")
    Polar.FFT(-1, "Space")
    W.FFT(-1, "Space")
    return W

def G_FirstOrder(Beta,G0, Sigma0, Sigma, map):
    G=weight.Weight("G.SmoothT", map, "TwoSpins", "AntiSymmetric")
    TauRange = range(G.Shape[TAU])
    SubRange=range(G.NSublat)
    SubList=[(a,b,c,d) for a in SubRange for b in SubRange for c in SubRange for d in SubRange]

    #make sure spin conservation on W0
    G0.FFT(1, "Space","Time")
    Sigma.FFT(1, "Space","Time")
    Sigma0.FFT(1, "Space")
    SpinList = [map.Spin2Index(UP,UP), map.Spin2Index(DOWN,DOWN)]
    for spG in SpinList:
        for e in SubList:
            subG0L = map.SublatIndex(e[0], e[1])
            subSigma0 = map.SublatIndex(e[1], e[2])
            subSigma = map.SublatIndex(e[1], e[2])
            subG0R = map.SublatIndex(e[2], e[3])
            subG = map.SublatIndex(e[0], e[3])
            G.Data[spG,subG,:,:]+=G0.Data[spG,subG0L,:,:] \
                    *Sigma.Data[spG,subSigma,:,:]*G0.Data[spG,subG0R,:,:]
            for tau in TauRange:
                G.Data[spG,subG,:,tau]+=G0.Data[spG,subG0L,:,tau] \
                        *Sigma0.Data[spG,subSigma,:]*G0.Data[spG,subG0R,:,tau]

    G.Data[:,:,:,:] *= (Beta/G.Shape[VOL]/(G.Shape[TAU]**2.0))

    G0.FFT(-1, "Space","Time")
    Sigma.FFT(-1, "Space","Time")
    Sigma0.FFT(-1, "Space")
    G.FFT(-1, "Space","Time")
    return G


def Sigma_FirstOrder(G, W, map):
    Sigma=weight.Weight("Sigma.SmoothT", map, "TwoSpins", "AntiSymmetric")

    TauRange = range(G.Shape[TAU])

    for spin1 in range(2):
        for spin2 in range(2):
            spinW = map.Spin4Index((spin1,spin2),(spin2,spin1))
            spinG = map.Spin2Index(spin2, spin2)
            spinSigma = map.Spin2Index(spin1, spin1)
            Sigma.Data[spinSigma, :, :, :]  \
                    = G.Data[spinG, :, :, :]\
                    *W.Data[spinW, :, :, :]

    return Sigma

def Sigma0_FirstOrder(G, W0, map):

    Sigma0=weight.Weight("Sigma.DeltaT", map, "TwoSpins", "AntiSymmetric")

    TauRange = range(G.Shape[TAU])
    for spin1 in range(2):
        for spin2 in range(2):
            spinW = map.Spin4Index((spin1,spin2),(spin2,spin1))
            spinG = map.Spin2Index(spin2, spin2)
            spinSigma = map.Spin2Index(spin1, spin1)
            for tau in TauRange:
                Sigma0.Data[spinSigma, :, :]  \
                        = G.Data[spinG, :, :, tau]\
                        *W0.Data[spinW, :, :]
    return Sigma0


def W_Dyson(Beta, W0,Polar,map):
    W=weight.Weight("W.SmoothT", map, "FourSpins", "Symmetric")
    W0.FFT(1, "Space")
    Polar.FFT(1, "Space", "Time")

    NSpin, NSub=W.NSpin, W.NSublat
    W0.Reshape("SPSUBSPSUB")
    W.Reshape("SPSUBSPSUB")
    Polar.Reshape("SPSUBSPSUB")

    JP=np.einsum("ijv,jkvt->ikvt",W0.Data, Polar.Data)
    #JP shape: NSpin*NSub,NSpin*NSub,Vol,Tau

    JP=(Beta/(W.Shape[TAU])**2.0) * JP
    for tau in range(W.Shape[TAU]):
        JP[:,:,:,tau] = JP[:,:,:,tau] * np.cos(tau*np.pi/W.Shape[TAU])

    I=np.eye(NSpin*NSub)
    W.Data=I[...,np.newaxis,np.newaxis]-JP
    W.Inverse();
    W.Data=np.einsum('ijvt,jkv->ikvt', W.Data,W0.Data)-W0.Data[...,np.newaxis]

    W.Reshape("SP2SUB2")
    W0.Reshape("SP2SUB2")
    Polar.Reshape("SP2SUB2")

    W.FFT(-1, "Space", "Time")
    Polar.FFT(-1, "Space", "Time")
    W0.FFT(-1, "Space")
    #print W.Data[map.Spin4Index((0,0),(0,0)),map.SublatIndex(0,0),0,:]
    return W

def G_Dyson(Beta, G0, Sigma0, Sigma, map):
    G=weight.Weight("G.SmoothT", map, "TwoSpins", "AntiSymmetric")
    G0.FFT(1, "Space", "Time")
    Sigma0.FFT(1, "Space")
    Sigma.FFT(1, "Space", "Time")

    NSpin, NSub=G.NSpin, G.NSublat

    G0.Reshape("SPSUBSPSUB")
    G.Reshape("SPSUBSPSUB")
    Sigma0.Reshape("SPSUBSPSUB")
    Sigma.Reshape("SPSUBSPSUB")

    G0Sigma0=np.einsum("ijvt,jkv->ikvt",G0.Data, Sigma0.Data)
    ####correction term
    for tau in range(G.Shape[TAU]):
        G0Sigma0[:,:,:,tau] = G0Sigma0[:,:,:,tau]*np.cos((tau+0.5)*np.pi/G.Shape[TAU])

    G0Sigma=np.einsum("ijvt,jkvt->ikvt",G0.Data, Sigma.Data)
    GS=(Beta/(G.Shape[VOL]*(G.Shape[TAU])**2.0)) * (G0Sigma0+G0Sigma)
    #GS shape: NSpin*NSub,NSpin*NSub,Vol,Tau

    I=np.eye(NSpin*NSub)
    G.Data=I[...,np.newaxis,np.newaxis]-GS
    G.Inverse();
    G.Data=np.einsum('ijvt,jkvt->ikvt', G.Data,G0.Data)

    G.Reshape("SP2SUB2")
    G0.Reshape("SP2SUB2")
    Sigma0.Reshape("SP2SUB2")
    Sigma.Reshape("SP2SUB2")

    G.FFT(-1, "Space", "Time")
    Sigma0.FFT(-1, "Space")
    Sigma.FFT(-1, "Space", "Time")
    G0.FFT(-1, "Space","Time")
    return G


