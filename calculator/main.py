#!usr/bin/env python
import numpy as np
import parameter as para
from weight import UP,DOWN,IN,OUT,TAU,SP,SUB,VOL
import weight

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
    W.SmoothT[:,:,:,:] = 0.0

def Calculate_Sigma_First_Order(G0, W, Sigma):
    Sigma.SmoothT[:,:,:,:] = 0.0

G0=weight.Weight("G", Beta, False)
G0.Load("../data/GW.npz")

W0=weight.Weight("W", Beta, True)
W0.Load("../data/GW.npz")

G=weight.Weight("G", Beta, False)
G.SmoothT = np.zeros(G0.SmoothT.shape, dtype=complex)

W=weight.Weight("W", Beta, True)
W.SmoothT = np.zeros(W0.DeltaT.shape+(G0.SmoothT.shape[-1],), dtype=complex)

Sigma=weight.Weight("Sigma", Beta, False)
Sigma.SmoothT = np.zeros(G.SmoothT.shape, dtype=complex)
Sigma.DeltaT = np.zeros(G.SmoothT.shape[0:-2], dtype=complex)

Polar=weight.Weight("Polar", Beta, True)
Polar.SmoothT = np.zeros(W.SmoothT.shape, dtype=complex)

Calculate_Polar_First_Order(G0, Polar)
map = G0.GetMap()
print G0.SmoothT[map.Spin2Index(UP,UP), map.SublatIndex(1,1),0,:]
print G0.SmoothT[map.Spin2Index(DOWN,DOWN), map.SublatIndex(1,1),0,::-1]
print Polar.SmoothT[map.Spin4Index((DOWN,UP),(UP,DOWN)), map.SublatIndex(1,1),0,:]

Calculate_W_First_Order(W0, Polar, W) 

Calculate_Sigma_First_Order(G0, W, Sigma)



