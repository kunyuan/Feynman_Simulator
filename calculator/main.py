#!usr/bin/env python
import numpy as np
import parameter as para
import calculator as calc
from weight import UP,DOWN,IN,OUT,TAU,SP,SUB,VOL
import weight
import model
from logger import *

para = para.Parameter()
para.Load("../data/infile/_in_DYSON_1")
Assert(para.Type=="DYSON", "The job type should be DYSON, not {0}".format(para.Type))
Beta = para.InitialBeta
Para={"NSublat":2, "L":para.L, "Beta": Beta, "MaxTauBin": 32}
map=weight.IndexMap(**Para)

if para.StartFromBare is True:
    Factory=model.BareFactory(map, para.Hopping, para.Interaction, 
                              para.ChemicalPotential, para.ExternalField)
    G0,W0=Factory.Build(para.Model, "Checkboard")
    #Factory.PlotModel()
else:
    G0=weight.Weight("G.SmoothT", map, "TwoSpins", "AntiSymmetric")
    G0.Load("../data/GW.npz")

W0=weight.Weight("W.DeltaT", map, "FourSpins", "Symmetric")
W0.Load("../data/GW.npz")

Polar=calc.Polar_FirstOrder(G0, map)
W=calc.W_Dyson(W0,Polar,map)

#W=calc.W_FirstOrder(W0, Polar,map) 
