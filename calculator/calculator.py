#!usr/bin/env python
import sys
import numpy as np
import parameter as para
from weight import UP,DOWN,IN,OUT,TAU,SP1,SUB1,SP2,SUB2,VOL,SPIN
import weight
from logger import *

def PlotTime(array, Beta):
    import matplotlib.pyplot as plt
    x=np.linspace(0, Beta, len(array))
    plt.plot(x,array,'-')
    plt.show()


def Sigma_FirstOrder(G, W, map):
    Sigma=weight.Weight("SmoothT", map, "TwoSpins", "AntiSymmetric")

    TauRange = range(map.MaxTauBin)

    for spin1 in range(2):
        for spin2 in range(2):
            spinW = (map.Spin2Index(spin1,spin2),map.Spin2Index(spin2,spin1))
            spinG = (spin2, spin2)
            spinSigma = (spin1, spin1)
            Sigma.Data[spinSigma[IN], :, spinSigma[OUT], :, :]  \
                    -= G.Data[spinG[IN], :, spinG[OUT], :, :]\
                    *W.Data[spinW[IN], :, spinW[OUT], :, :]
    return Sigma

def Sigma0_FirstOrder(G, W0, map):
    Sigma0=weight.Weight("DeltaT", map, "TwoSpins", "AntiSymmetric")

    for spin1 in range(2):
        for spin2 in range(2):
            #############G(tau==-0)
            spinW = (map.Spin2Index(spin1,spin2),map.Spin2Index(spin2,spin1))
            spinG = (spin2, spin2)
            spinSigma = (spin1, spin1)
            Sigma0.Data[spinSigma[IN], :, spinSigma[OUT], :, :]  \
                    -= G.Data[spinG[IN], :, spinG[OUT], :, :, -1]\
                    *W0.Data[spinW[IN], :, spinW[OUT], :, :]
    return Sigma0


def Polar_FirstOrder(G, map):
    Polar=weight.Weight("SmoothT", map, "FourSpins","Symmetric")
    NSublat = G.NSublat
    SubList=[(a,b) for a in range(NSublat) for b in range(NSublat)]
    SpList=[(a,b) for a in range(SPIN) for b in range(SPIN)]
    for spin1,spin2 in SpList:
        spinPolar = ((spin1,spin2),(spin2,spin1))
        spinG1 = (spin1, spin1)
        spinG2 = (spin2, spin2)
        for subA,subB in SubList:
            Polar.Data[map.Spin2Index(*spinPolar[IN]),subA, \
                    map.Spin2Index(*spinPolar[OUT]),subB,:,:]\
                    = (-1.0)*G.Data[spinG1[IN], subB, spinG1[OUT], subA, :, ::-1]  \
                    *G.Data[spinG2[IN], subA, spinG2[OUT], subB, :, :]
    return Polar

def W_FirstOrder(Beta,W0, Polar, map):
    W=weight.Weight("SmoothT", map, "FourSpins", "Symmetric")
    TauRange = range(map.MaxTauBin)
    SubRange=range(W.NSublat)
    SubList=[(a,b,c,d) for a in SubRange for b in SubRange for c in SubRange for d in SubRange]
    SpinList=[(Wtuple,Polartuple) for Wtuple in map.GetConservedSpinTuple("FourSpins") \
                                  for Polartuple in map.GetConservedSpinTuple("FourSpins") \
                                  if map.IsConserved(4, (Wtuple[IN], Polartuple[IN]))] 
    #make sure spin conservation on W0

    W0.FFT(1, "Space")
    Polar.FFT(1, "Space")

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
    
    W0.FFT(-1, "Space")
    Polar.FFT(-1, "Space")
    W.FFT(-1, "Space")
    return W

def W_Dyson(Beta, W0, Polar, map):
    W=weight.Weight("SmoothT", map, "FourSpins", "Symmetric")
    W0.FFT(1, "Space")
    Polar.FFT(1, "Space", "Time")

    NSpin, NSub=W.NSpin, W.NSublat

    JP=np.einsum("ijklv,klmnvt->ijmnvt",W0.Data, Polar.Data)
    #JP shape: NSpin,NSub,NSpin,NSub,Vol,Tau

    JP *= Beta/map.MaxTauBin
    for tau in range(map.MaxTauBin):
        JP[:,:,:,:,:,tau] = JP[:,:,:,:,:,tau] * np.cos(tau*np.pi/map.MaxTauBin)

    I=np.eye(NSpin*NSub).reshape([NSpin,NSub,NSpin,NSub])
    W.Data=I[...,np.newaxis,np.newaxis]-JP
    W.Inverse();
    W.Data = np.einsum('ijklvt,klmnv->ijmnvt', W.Data,W0.Data)
    W.Data = map.MaxTauBin/Beta*(W.Data - W0.Data[...,np.newaxis])

    W.FFT(-1, "Space", "Time")
    Polar.FFT(-1, "Space", "Time")
    W0.FFT(-1, "Space")
    return W

def G_Dyson(Beta, G0, Sigma0, Sigma, map):
    G=weight.Weight("SmoothT", map, "TwoSpins", "AntiSymmetric")
    G0.FFT(1, "Space", "Time")
    Sigma0.FFT(1, "Space")
    Sigma.FFT(1, "Space", "Time")

    NSpin, NSub=G.NSpin, G.NSublat

    G0Sigma0=np.einsum("ijklvt,klmnv->ijmnvt",G0.Data, Sigma0.Data)
    G0Sigma=np.einsum("ijklvt,klmnvt->ijmnvt",G0.Data, Sigma.Data)

    ####correction term
    for tau in range(map.MaxTauBin):
        G0Sigma0[:,:,:,:,:,tau] = G0Sigma0[:,:,:,:,:,tau]*np.cos(np.pi*map.IndexToTau(tau)/Beta)

    GS  = Beta/map.MaxTauBin*(Beta/map.MaxTauBin*G0Sigma + G0Sigma0)
    #GS shape: NSpin,NSub,NSpin,NSub,Vol,Tau

    I=np.eye(NSpin*NSub).reshape([NSpin,NSub,NSpin,NSub])
    G.Data=I[...,np.newaxis,np.newaxis]-GS
    G.Inverse();
    G.Data=np.einsum('ijklvt,klmnvt->ijmnvt', G.Data,G0.Data)

    G.FFT(-1, "Space", "Time")
    Sigma0.FFT(-1, "Space")
    Sigma.FFT(-1, "Space", "Time")
    G0.FFT(-1, "Space","Time")
    return G


def Calculate_Chi(Beta, W0, Polar, map):

    Chi=weight.Weight("SmoothT", map, "FourSpins", "Symmetric")

    W0.FFT(1, "Space")
    Polar.FFT(1, "Space", "Time")

    NSpin, NSub=Chi.NSpin, Chi.NSublat

    JP=np.einsum("ijklv,klmnvt->ijmnvt",W0.Data, Polar.Data)
    #JP shape: NSpin,NSub,NSpin,NSub,Vol,Tau

    JP *= Beta/map.MaxTauBin
    for tau in range(map.MaxTauBin):
        JP[:,:,:,:,:,tau] = JP[:,:,:,:,:,tau] * np.cos(tau*np.pi/map.MaxTauBin)

    I=np.eye(NSpin*NSub).reshape([NSpin,NSub,NSpin,NSub])
    Chi.Data=I[...,np.newaxis,np.newaxis]-JP
    Chi.Inverse();
    Chi.Data = np.einsum('ijklvt,klmnvt->ijmnvt', Polar.Data, Chi.Data)
    Chi.Data *=-1.5

    Chi.FFT(-1, "Space", "Time")
    Polar.FFT(-1, "Space", "Time")
    W0.FFT(-1, "Space")
    return Chi
