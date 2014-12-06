#!usr/bin/env python
import numpy as np
import parameter as para
from weight import UP,DOWN,IN,OUT,TAU,SP,SUB,VOL
import weight
from logger import *

para = para.Parameter()
para.Load("../data/infile/_in_DYSON_1")
Beta = para.InitialBeta


def Calculate_Polar_First_Order(G, Polar):
    map = G.GetMap()
    NSublat = 2   ######
    Ntau = G.SmoothT.shape[TAU]
    for spin1 in range(2):
        for spin2 in range(2):
            spinPolar = map.Spin4Index((spin1,spin2),(spin2,spin1))
            spinG1 = map.Spin2Index(spin1, spin1)
            spinG2 = map.Spin2Index(spin2, spin2)
            for subA in range(NSublat):
                for subB in range(NSublat):
                    subA2B = map.SublatIndex(subA, subB)
                    subB2A = map.SublatIndex(subB, subA)
                    Polar.SmoothT[spinPolar, subA2B, :, :]  \
                            = (-1.0)*G.SmoothT[spinG1, subB2A, :, ::-1]\
                            *G.SmoothT[spinG2, subA2B, :, :]


def Calculate_W_First_Order(W0, Polar, W):
    map = W0.GetMap()
    NSublat = 2 ########
    Ntau = Polar.SmoothT.shape[TAU]
    Nvol = Polar.SmoothT.shape[VOL]
    rsub=range(NSublat)
    subList=[(a,b,c,d) for a in rsub for b in rsub for c in rsub for d in rsub]
    a=W0.DeltaT.copy()
    b=Polar.SmoothT.copy()
    Assert(np.allclose(b[map.Spin4Index((UP,UP),(UP,UP)),...], b[map.Spin4Index((DOWN,DOWN),(DOWN,DOWN)),...]), "test polar")
    
    W0.fftSpace(1)
    Polar.fftSpace(1)
    Assert(np.allclose(Polar.SmoothT[map.Spin4Index((UP,UP),(UP,UP)),...], Polar.SmoothT[map.Spin4Index((DOWN,DOWN),(DOWN,DOWN)),...]), "test polar")
    
    for spWt in map.GetLegalSpinTuple(4):
        for spPolart in map.GetLegalSpinTuple(4):
            if (not map.IsLegal(4, (spWt[IN],spPolart[IN]))):
                continue
            spW0L=map.Spin4Index(spWt[IN], spPolart[IN])
            spW0R=map.Spin4Index(spPolart[OUT], spWt[IN])
            spW = map.Spin4Index(*spWt)
            spPolar = map.Spin4Index(*spPolart)
            for e in subList:
                subW0L = map.SublatIndex(e[0], e[1])
                subPolar = map.SublatIndex(e[1], e[2])
                subW0R = map.SublatIndex(e[2], e[3])
                subW = map.SublatIndex(e[0], e[3])
                for tau in range(Ntau):
                    W.SmoothT[spW,subW,:,tau]+=W0.DeltaT[spW0L,subW0L,:] \
                        *Polar.SmoothT[spPolar,subPolar,:,tau] \
                        *W0.DeltaT[spW0R,subW0R,:]
    W0.fftSpace(-1)
    Polar.fftSpace(-1)
    W.fftSpace(-1)

def Calculate_Sigma_First_Order(G0, W, Sigma):
    Sigma.SmoothT[:,:,:,:] = 0.0

G0=weight.Weight("G", Beta, para.L, False)
G0.Load("../data/GW.npz")

W0=weight.Weight("W", Beta, para.L, True)
W0.Load("../data/GW.npz")

G=weight.Weight("G", Beta, para.L, False)
G.SmoothT = np.zeros(G0.SmoothT.shape, dtype=complex)

W=weight.Weight("W", Beta, para.L, True)
W.SmoothT = np.zeros(W0.DeltaT.shape+(G0.SmoothT.shape[-1],), dtype=complex)

Sigma=weight.Weight("Sigma", Beta, para.L, False)
Sigma.SmoothT = np.zeros(G.SmoothT.shape, dtype=complex)
Sigma.DeltaT = np.zeros(G.SmoothT.shape[0:-2], dtype=complex)

Polar=weight.Weight("Polar", Beta, para.L, True)
Polar.SmoothT = np.zeros(W.SmoothT.shape, dtype=complex)

Calculate_Polar_First_Order(G0, Polar)

map = G0.GetMap()
print G0.SmoothT[map.Spin2Index(UP,UP), map.SublatIndex(1,1),0,:]
print G0.SmoothT[map.Spin2Index(DOWN,DOWN), map.SublatIndex(1,1),0,::-1]
print Polar.SmoothT[map.Spin4Index((DOWN,UP),(UP,DOWN)), map.SublatIndex(1,1),0,:]
print Polar.SmoothT[map.Spin4Index((UP,UP),(UP,UP)), map.SublatIndex(1,1),0,:]

Calculate_W_First_Order(W0, Polar, W) 
print W.SmoothT[map.Spin4Index((DOWN,DOWN),(DOWN,DOWN)), 0,0,:]
print W.SmoothT[map.Spin4Index((UP,UP),(UP,UP)), 0,0,:]

Calculate_Sigma_First_Order(G0, W, Sigma)

