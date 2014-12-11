#!usr/bin/env python
import numpy as np
import parameter as para
import calculator as calc
from weight import UP,DOWN,IN,OUT,TAU,SP,SUB,VOL,OneSpin,TwoSpin,Symmetric,AntiSymmetric
import weight
from logger import *

para = para.Parameter()
para.Load("../data/infile/_in_DYSON_1")
Beta = para.InitialBeta
DeltaTPara={"NSublat":2, "L":para.L}
SmoothTPara=DeltaTPara.copy()
SmoothTPara.update({"Beta": Beta, "MaxTauBin": 32})
map=weight.IndexMap(**SmoothTPara)

G0=weight.Weight("G", OneSpin, SmoothTPara, AntiSymmetric)
G0.Load("../data/GW.npz")

W0=weight.Weight("W", TwoSpin, DeltaTPara, Symmetric)
W0.Load("../data/GW.npz")

Polar=calc.Polar_FirstOrder(G0, map)

#print G0.Data[map.Spin2Index(UP,UP), map.SublatIndex(1,1),0,:]
#print G0.Data[map.Spin2Index(DOWN,DOWN), map.SublatIndex(1,1),0,::-1]
#print Polar.Data[map.Spin4Index((DOWN,UP),(UP,DOWN)), map.SublatIndex(1,1),0,:]
#print Polar.Data[map.Spin4Index((UP,UP),(UP,UP)), map.SublatIndex(1,1),0,:]

W=calc.W_FirstOrder(W0, Polar,map) 
#print W.Data[map.Spin4Index((DOWN,DOWN),(DOWN,DOWN)), 0,0,:]
#print W.Data[map.Spin4Index((UP,UP),(UP,UP)), 0,0,:]

W=calc.W_Dyson(W0,Polar,map)
#Sigma=calc.Sigma_FirstOrder(G0, W, weight.AntiSymmetric)

