#!/usr/bin/env python
import numpy as np
import calculator as calc
import lattice as lat
import collect
from weight import UP,DOWN,IN,OUT
from logger import *
import os, sys, model, weight, parameter, plot, argparse

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--PID", help="use PID to find the input file")
parser.add_argument("-f", "--file", help="use file path to find the input file")
args = parser.parse_args()
if args.PID:
    InputFile=os.path.join(workspace, "infile/_in_DYSON_"+str(args.PID))
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
print "G0=\n", G0.Data[UP,0,UP,0,0,:]
#Factory.Plot()

if para["Job"]["StartFromBare"] is True or os.path.exists(WeightFile+".pkl") is False:
    #start from bare
    log.info("Start from G0 and W0 to do dyson...")
    G=G0.Copy()
    for sub in range(2):
        print G.Data[UP,sub,UP,sub,0,-1]-G.Data[DOWN, sub, DOWN, sub, 0, -1]
    W=weight.Weight("SmoothT", Map, "FourSpins", "Symmetric")
    #ChiTensor = calc.Calculate_ChiTensor(W0, Polar, Map)
    #Chi, _=calc.Calculate_Chi(ChiTensor, Map)
    #plot.PlotSpatial(Chi, Lat, 0, 0)
    #sys.exit(0)

    for i in range(10):
        log.info("Round {0}...".format(i))
        Sigma=calc.Sigma_FirstOrder(G, W, Map)
        Sigma0=calc.Sigma0_FirstOrder(G, W0, Map)
        Polar=calc.Polar_FirstOrder(G, Map)

        #######DYSON FOR W AND G###########################
        G = calc.G_Dyson(G0, Sigma0, Sigma, Map)
        W = calc.W_Dyson(W0, Polar, Map)
        ###################################################
        calc.Check_Denorminator(W0, Polar, Map)

else:
    #########READ G,SIGMA,POLAR; CALCULATE SIGMA0 #################
    log.info("Load Sigma/Polar and G to do dyson...")
    data=IO.LoadBigDict(WeightFile)
    G=weight.Weight("SmoothT", Map, "TwoSpins", "AntiSymmetric").FromDict(data["G"])
    MaxOrder=para["Dyson"]["Order"]
    SigmaMC, PolarMC=collect.CollectStatis(Map, MaxOrder)
    Sigma,Polar=collect.UpdateWeight(SigmaMC, PolarMC, paraDyson["ErrorThreshold"], paraDyson["OrderAccepted"])

    Sigma0=calc.Sigma0_FirstOrder(G, W0, Map)
    #######DYSON FOR W AND G###########################
    G = calc.G_Dyson(G0, Sigma0, Sigma, Map)
    W = calc.W_Dyson(W0, Polar,Map)
    ###################################################

######Calculate Chi ###############################
ChiTensor = calc.Calculate_ChiTensor(W0, Polar, Map)
Chi, _=calc.Calculate_Chi(ChiTensor, Map)

##########OUTPUT AND FILE SAVE ####################
spinUP=Map.Spin2Index(UP,UP)
print "Polar=\n", Polar.Data[spinUP,0,spinUP,0,0,:]
print "W=\n", W.Data[spinUP,0,spinUP,0,0,:]
print "G=\n", G.Data[UP,0,UP,0,0,:]
print "Sigma=\n", Sigma.Data[UP,0,UP,0,0,:]
print "Polar=\n", Polar.Data[spinUP,0,spinUP,0,0,:]
print "Chi=\n", Chi.Data[0,0,0,0,0,:]

data={}
data["G"]=G.ToDict()
data["W"]=W.ToDict()
data["W"].update(W0.ToDict())
data["Sigma"]=Sigma.ToDict()    ####ForTest
data["Polar"]=Polar.ToDict()    ####ForTest
data["Chi"]=Chi.ToDict()        ####ForTest

IO.SaveBigDict(WeightFile, data)
###################################################
plot.PlotSpatial(Chi, Lat, 0, 0)
