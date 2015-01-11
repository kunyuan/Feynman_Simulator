#!/usr/bin/env python
import numpy as np
import calculator as calc
import lattice as lat
import collect
from weight import UP,DOWN,IN,OUT,TAU,SP1,SUB1,SP2,SUB2,VOL
from logger import *
import os, sys, model, weight, parameter, plot, argparse

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--PID", help="use PID to find the input file")
parser.add_argument("-f", "--file", help="use file path to find the input file")
args = parser.parse_args()
if args.PID:
    InputFile=workspace+"infile/_in_DYSON_"+str(args.PID)
elif args.file:
    InputFile=os.path.abspath(args.file)
else:
    Assert(False, "Do not understand the argument!")

para=parameter.Load(InputFile)
WeightFile=para["Job"]["WeightFile"]
WeightPara={"NSublat": para["Lattice"]["NSublat"], "L":para["Lattice"]["L"],
            "Beta": para["Tau"]["Beta"], "MaxTauBin": para["Tau"]["MaxTauBin"]}
Map=weight.IndexMap(**WeightPara)
Lat=lat.Lattice(para["Lattice"]["Name"], Map)

##########INITIALIZATION ##########################
Factory=model.BareFactory(Map, para["Model"])
G0,W0=Factory.Build(para["Model"]["Name"], para["Lattice"]["Name"])
#Factory.Plot()

if para["Job"]["StartFromBare"] is True or os.path.exists(WeightFile+".pkl") is False:
    #start from bare
    log.info("Start from G0 and W0 to do dyson...")
    G=G0.Copy()
    W=weight.Weight("SmoothT", Map, "FourSpins", "Symmetric")
    for i in range(10):
        log.info("Round {0}...".format(i))
        Polar=calc.Polar_FirstOrder(G, Map)
        Sigma=calc.Sigma_FirstOrder(G, W, Map)
        Sigma0=calc.Sigma0_FirstOrder(G, W0, Map)
        #######DYSON FOR W AND G###########################
        W = calc.W_Dyson(W0, Polar,map)
        G = calc.G_Dyson(G0, Sigma0, Sigma, Map)
        ###################################################
else:
    #########READ G,SIGMA,POLAR; CALCULATE SIGMA0 #################
    log.info("Load Sigma/Polar and G to do dyson...")
    data=IO.LoadBigDict(WeightFile)
    G=weight.Weight("SmoothT", Map, "TwoSpins", "AntiSymmetric").FromDict(data["G"])
    #Sigma=weight.Weight("SmoothT", map, "TwoSpins", "AntiSymmetric").FromDict(data["Sigma"])
    #Polar=weight.Weight("SmoothT", map, "FourSpins", "Symmetric").FromDict(data["Polar"])
    paraDyson=para["Dyson"]
    Sigma, Polar=collect.CollectStatis(Map, paraDyson["Order"],
                                       paraDyson["ErrorThreshold"], paraDyson["OrderAccepted"])

    Sigma0=calc.Sigma0_FirstOrder(G, W0, Map)
    #######DYSON FOR W AND G###########################
    W = calc.W_Dyson(W0, Polar,Map)
    G = calc.G_Dyson(G0, Sigma0, Sigma, Map)
    ###################################################

######Calculate Chi ###############################
Chi = calc.Calculate_Chi(W0, Polar, Map)

##########OUTPUT AND FILE SAVE ####################
spinUP=Map.Spin2Index(UP,UP)
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


