#!/usr/bin/env python
import numpy as np
import calculator as calc
import lattice as lat
from weight import UP,DOWN,IN,OUT,TAU,SP1,SUB1,SP2,SUB2,VOL
from logger import *
import os, sys, model, weight, parameter, plot

para=parameter.Load(os.path.abspath(sys.argv[1]))
WeightFile=para["Job"]["WeightFile"]
WeightPara={"NSublat": para["Lattice"]["NSublat"], "L":para["Lattice"]["L"],
            "Beta": para["Tau"]["Beta"], "MaxTauBin": para["Tau"]["MaxTauBin"]}
map=weight.IndexMap(**WeightPara)
Lat=lat.Lattice(para["Lattice"]["Name"], map)

##########INITIALIZATION ##########################
Factory=model.BareFactory(map, para["Model"])
G0,W0=Factory.Build(para["Model"]["Name"], para["Lattice"]["Name"])
#Factory.Plot()

if para["Job"]["DoesLoad"] is False:
    G=G0.Copy()
    W=weight.Weight("SmoothT", map, "FourSpins", "Symmetric")
    for i in range(10):
        log.info("Round {0}...".format(i))
        Polar=calc.Polar_FirstOrder(G, map)
        Sigma=calc.Sigma_FirstOrder(G, W, map)
        Sigma0=calc.Sigma0_FirstOrder(G, W0, map)
        #######DYSON FOR W AND G###########################
        W = calc.W_Dyson(W0, Polar,map)
        G = calc.G_Dyson(G0, Sigma0, Sigma, map)
        ###################################################

else:
    #########READ G,SIGMA,POLAR; CALCULATE SIGMA0 #################
    data=IO.LoadBigDict(WeightFile)
    G=weight.Weight("SmoothT", map, "TwoSpins", "AntiSymmetric").FromDict(data["G"])
    Sigma=weight.Weight("SmoothT", map, "TwoSpins", "AntiSymmetric").FromDict(data["Sigma"])
    Polar=weight.Weight("SmoothT", map, "FourSpins", "Symmetric").FromDict(data["Polar"])

    Sigma0=calc.Sigma0_FirstOrder(G, W0, map)
    #######DYSON FOR W AND G###########################
    W = calc.W_Dyson(W0, Polar,map)
    G = calc.G_Dyson(G0, Sigma0, Sigma, map)
    ###################################################

######Calculate Chi ###############################
Chi = calc.Calculate_Chi(W0, Polar, map)

##########OUTPUT AND FILE SAVE ####################
spinUP=map.Spin2Index(UP,UP)
#print "Polar=\n", Polar.Data[spinUP,0,spinUP,0,0,:]
print "W=\n", W.Data[spinUP,0,spinUP,0,0,:]
print "G=\n", G.Data[UP,0,UP,0,0,:]
print "Chi=\n", Chi.Data[spinUP,0,spinUP,0,0,:]

print WeightFile
data={}
data["G"]=G.ToDict()
data["W"]=W.ToDict()
data["W"].update(W0.ToDict())

data["Sigma"]=Sigma.ToDict()    ####ForTest
data["Polar"]=Polar.ToDict()    ####ForTest
data["Chi"]=Chi.ToDict()        ####ForTest

IO.SaveBigDict(WeightFile, data)
###################################################
plot.PlotSpatial(Chi, Lat, spinUP, spinUP)


