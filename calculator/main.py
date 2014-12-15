#!usr/bin/env python
import numpy as np
import parameter as para
import calculator as calc
from weight import UP,DOWN,IN,OUT,TAU,SP,SUB,VOL
import weight
import model
import lattice as lat
from logger import *

prefix="../data/"
para = para.Parameter()
para.Load("../data/infile/_in_DYSON_1")
Assert(para.Type=="DYSON", "The job type should be DYSON, not {0}".format(para.Type))
Beta = para.InitialBeta
Lat=lat.Lattice(para.Lattice, para.L)
WeightPara={"NSublat":Lat.NSublat, "L":para.L, "Beta": Beta, "MaxTauBin": 32}
map=weight.IndexMap(**WeightPara)

if para.StartFromBare is True:
    Factory=model.BareFactory(map, para.Hopping, para.Interaction, 
                              para.ChemicalPotential, para.ExternalField)
    G0,W0=Factory.Build(para.Model, para.Lattice)
    G0.Save(prefix+para.WeightFile,"w")
    W0.Save(prefix+para.WeightFile,"a")
    Factory.Plot()
else:
    G0=weight.Weight("G.SmoothT", map, "TwoSpins", "AntiSymmetric")
    G0.Load("../data/GW.npz")
    W0=weight.Weight("W.DeltaT", map, "FourSpins", "Symmetric")
    W0.Load("../data/GW.npz")

Polar=calc.Polar_FirstOrder(G0, map)
print "Polar", Polar.Data[map.Spin4Index((DOWN,UP),(UP,DOWN)), 0,0,:]

W=calc.W_FirstOrder(Beta, W0, Polar,map) 
print W.Data[map.Spin4Index((DOWN,DOWN),(DOWN,DOWN)), 0,0,:]

Sigma=calc.Sigma_FirstOrder(G0, W, map)
Sigma0=calc.Sigma0_FirstOrder(G0, W0, map)

W = calc.W_Dyson(Beta, W0,Polar,map)
print W.Data[map.Spin4Index((DOWN,DOWN),(DOWN,DOWN)), 0,0,:]

G = calc.G_FirstOrder(Beta,G0, Sigma0, Sigma, map) 
print G0.Data[map.Spin2Index(UP,UP),0,0,:]+G.Data[map.Spin2Index(UP,UP), 0,0,:]

G = calc.G_Dyson(Beta, G0, Sigma0, Sigma, map)
print G.Data[map.Spin2Index(UP,UP), 0,0,:]

#W=calc.W_FirstOrder(W0, Polar,map) 
