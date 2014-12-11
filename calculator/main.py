#!usr/bin/env python
import numpy as np
import parameter as para
import calculator as calc
from weight import UP,DOWN,IN,OUT,TAU,SP,SUB,VOL
import weight
from logger import *

para = para.Parameter()
para.Load("../data/infile/_in_DYSON_1")
Beta = para.InitialBeta
Para={"NSublat":2, "L":para.L, "Beta": Beta, "MaxTauBin": 32}
map=weight.IndexMap(**Para)

G0=weight.Weight("G.SmoothT", map, "OneSpin", "AntiSymmetric")
G0.Load("../data/GW.npz")

W0=weight.Weight("W.DeltaT", map, "TwoSpins", "Symmetric")
W0.Load("../data/GW.npz")

Polar=calc.Polar_FirstOrder(G0, map)

#print G0.Data[map.Spin2Index(UP,UP), map.SublatIndex(1,1),0,:]
#print G0.Data[map.Spin2Index(DOWN,DOWN), map.SublatIndex(1,1),0,::-1]
#print Polar.Data[map.Spin4Index((DOWN,UP),(UP,DOWN)), map.SublatIndex(1,1),0,:]
#print Polar.Data[map.Spin4Index((UP,UP),(UP,UP)), map.SublatIndex(1,1),0,:]

W=calc.W_FirstOrder(W0, Polar,map) 
#print W.Data[map.Spin4Index((DOWN,DOWN),(DOWN,DOWN)), 0,0,:]
#print W.Data[map.Spin4Index((UP,UP),(UP,UP)), 0,0,:]

W=calc.W_FirstOrder(W0, Polar) 
print W.SmoothT[map.Spin4Index((DOWN,DOWN),(DOWN,DOWN)), 0,0,:]
print W.SmoothT[map.Spin4Index((UP,UP),(UP,UP)), 0,0,:]

Sigma=calc.Sigma_FirstOrder(G0, W0, W, weight.AntiSymmetric)
print Sigma.SmoothT[map.Spin2Index(UP,UP), map.SublatIndex(1,1),0,:]
print Sigma.DeltaT[map.Spin2Index(UP,UP), map.SublatIndex(1,1),0]

W=calc.W_Dyson(W0,Polar,map)

