#!usr/bin/env python
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
    Polar=weight.Weight("Polar.SmoothT", map, "TwoSpins","Symmetric")
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
    W=weight.Weight("W.SmoothT", map, "TwoSpins", "Symmetric")
    W0.FFT(1, "Space")
    Polar.FFT(1, "Space", "Time")

    NSpin, NSub=W.NSpin, W.NSublat
    W0.Reshape("SPSUBSPSUB")
    W.Reshape("SPSUBSPSUB")
    Polar.Reshape("SPSUBSPSUB")
    JP=np.einsum("ijv,jkvt->ikvt",W0.Data, Polar.Data)
    #JP shape: NSpin*NSub,NSpin*NSub,Vol,Tau
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

def Sigma_FirstOrder(G0, W, map):
    Sigma=weight.Weight("Sigma.SmoothT", map, "OneSpins", "AntiSymmetric")
    return Sigma

def W_FirstOrder(W0, Polar, map):
    W=weight.Weight("W.SmoothT", map, "TwoSpins", "Symmetric")
    TauRange = range(W.Shape[TAU])
    SubRange=range(W.NSublat)
    SubList=[(a,b,c,d) for a in SubRange for b in SubRange for c in SubRange for d in SubRange]
    SpinList=[(Wtuple,Polartuple) for Wtuple in map.GetConservedSpinTuple(4) \
                                  for Polartuple in map.GetConservedSpinTuple(4) \
                                  if map.IsConserved(4, (Wtuple[IN], Polartuple[IN]))] #make sure spin conservation on W0
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
    W0.FFT(-1, "Space")
    Polar.FFT(-1, "Space")
    W.FFT(-1, "Space")
    return W
