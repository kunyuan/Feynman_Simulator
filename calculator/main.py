#!usr/bin/env python
import numpy as np
import parameter as para
import calculator as calc
from weight import UP,DOWN,IN,OUT,TAU,SP1,SUB1,SP2,SUB2,VOL
import weight
import model
import lattice as lat
from logger import *
import sys
sys.path.append("../") #add the root dir into PYTHONPATH
import IO

prefix="../data/"
para = para.Parameter()
para.Load("../data/infile/_in_DYSON_1")
Assert(para.Type=="DYSON", "The job type should be DYSON, not {0}".format(para.Type))
Beta = para.InitialBeta
WeightPara={"NSublat":para.NSublat, "L":para.L, "Beta": Beta, "MaxTauBin": para.MaxTauBin}
map=weight.IndexMap(**WeightPara)

##########INITIALIZATION ##########################
Factory=model.BareFactory(map, para.Hopping, para.Interaction, 
                            para.ChemicalPotential, para.ExternalField)
G0,W0=Factory.Build(para.Model, para.Lattice)
data={}
data["G"]=G0.ToDict()
data["W"]=W0.ToDict()
IO.SaveBigDict(prefix+para.WeightFile, data)
#Factory.Plot()

if para.StartFromBare is True:
    #####FIRST ORDER CALCULATION FOR POLAR, W AND SIGMA############
    Polar=calc.Polar_FirstOrder(G0, map)
    W = calc.W_Dyson(Beta, W0, Polar,map)

    Sigma=calc.Sigma_FirstOrder(G0, W, map)
    Sigma0=calc.Sigma0_FirstOrder(G0, W0, map)

else:
    #########READ G,SIGMA,POLAR; CALCULATE SIGMA0 #################
    data=IO.LoadBigDict(prefix+para.WeightFile)
    G=weight.Weight("SmoothT", map, "TwoSpins", "AntiSymmetric")
    #G.Load(prefix+para.WeightFile)
    G.FromDict(data["G"])
    Sigma0=calc.Sigma0_FirstOrder(G, W0, map)

    Sigma=weight.Weight("SmoothT", map, "TwoSpins", "AntiSymmetric")
    #Sigma.Load(prefix+para.WeightFile)
    Sigma.FromDict(data["Sigma"])

    Polar=weight.Weight("SmoothT", map, "FourSpins", "Symmetric")
    #Polar.Load(prefix+para.WeightFile)
    Polar.FromDict(data["Polar"])

#######DYSON FOR W AND G###########################
W = calc.W_Dyson(Beta, W0, Polar,map)
G = calc.G_Dyson(Beta, G0, Sigma0, Sigma, map)
###################################################


##########OUTPUT AND FILE SAVE ####################
spinUP=map.Spin2Index(UP,UP)
print W.Data[spinUP,0,spinUP,0,0,:]
print G.Data[UP,0,UP,0,0,:]

data={}
data["G"]=G.ToDict()
data["W"]=W.ToDict()
IO.SaveBigDict(prefix+para.WeightFile, data)


###################################################

