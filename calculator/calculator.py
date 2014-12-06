#!usr/bin/env python
import numpy as np
import parameter as para
from weight import UP,DOWN,IN,OUT,TAU,SP,SUB,VOL
import weight
from logger import *

def Polar_FirstOrder(G):
    map = G.GetMap()
    Polar=weight.Weight("Polar", G.Beta, G.L, weight.Symmetric)
    Polar.SetShape(map.GetShape(4))
    Polar.SetZeros("SmoothT")
    NSublat = G.NSublattice
    SubList=[(map.SublatIndex(a,b),map.SublatIndex(b,a)) for a in range(NSublat) for b in range(NSublat)]
    for spin1 in range(2):
        for spin2 in range(2):
            spinPolar = map.Spin4Index((spin1,spin2),(spin2,spin1))
            spinG1 = map.Spin2Index(spin1, spin1)
            spinG2 = map.Spin2Index(spin2, spin2)
            for subA2B,subB2A in SubList:
                Polar.SmoothT[spinPolar, subA2B, :, :]  \
                        = (-1.0)*G.SmoothT[spinG1, subB2A, :, ::-1]\
                        *G.SmoothT[spinG2, subA2B, :, :]
    return Polar

def W_FirstOrder(W0, Polar):
    map = W0.GetMap()
    W=weight.Weight("W", W0.Beta, W0.L, weight.Symmetric)
    W.SetShape(map.GetShape(4))
    W.SetZeros("SmoothT")
    TauRange = range(W0.TauBinMax)
    SubRange=range(W0.NSublattice)
    SubList=[(a,b,c,d) for a in SubRange for b in SubRange for c in SubRange for d in SubRange]
    SpinList=[(Wtuple,Polartuple) for Wtuple in map.GetLegalSpinTuple(4) \
                                  for Polartuple in map.GetLegalSpinTuple(4) \
                                  if map.IsLegal(4, (Wtuple[IN], Polartuple[IN]))] #make sure spin conservation on W0
    print "\n".join([str(e) for e in SpinList])

    W0.fftSpace(1)
    Polar.fftSpace(1)
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
                W.SmoothT[spW,subW,:,tau]+=W0.DeltaT[spW0L,subW0L,:] \
                    *Polar.SmoothT[spPolar,subPolar,:,tau]*W0.DeltaT[spW0R,subW0R,:]
    W0.fftSpace(-1)
    Polar.fftSpace(-1)
    W.fftSpace(-1)
    return W

def Sigma_FirstOrder(G0, W, IsSymmetric=weight.AntiSymmetric):
    map = G0.GetMap()
    Sigma=weight.Weight("Sigma", G0.Beta, G0.L, IsSymmetric)
    Sigma.SetShape(map.GetShape(2))
    Sigma.SetZeros("SmoothT")
    return Sigma
