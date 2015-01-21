#!/usr/bin/env python
import numpy as np
import calculator as calc
import lattice as lat
import collect
from weight import UP,DOWN,IN,OUT
from logger import *
import os, sys, model, weight, parameter, plot, argparse, time, traceback

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--PID", help="use PID to find the input file")
parser.add_argument("-f", "--file", help="use file path to find the input file")
parser.add_argument("-c", "--collect", action="store_true", help="collect all the _statis.pkl file into statis_total.pkl")
args = parser.parse_args()
if args.PID:
    InputFile=os.path.join(workspace, "infile/_in_DYSON_"+str(args.PID))
elif args.file:
    InputFile=os.path.abspath(args.file)
else:
    Assert(False, "Do not understand the argument!")

job, para=parameter.Load(InputFile)
WeightFile=job["WeightFile"]
WeightPara={"NSublat": para["Lattice"]["NSublat"], "L":para["Lattice"]["L"],
            "Beta": para["Tau"]["Beta"], "MaxTauBin": para["Tau"]["MaxTauBin"]}
Map=weight.IndexMap(**WeightPara)
Lat=lat.Lattice(para["Lattice"]["Name"], Map)


if args.collect:
    MaxOrder=para["Dyson"]["Order"]
    SigmaMC, PolarMC=collect.CollectStatis(Map, MaxOrder)
    data ={}
    data["Sigma"] = {"Histogram": SigmaMC.ToDict()}
    data["Polar"] = {"Histogram": PolarMC.ToDict()}
    IO.SaveBigDict(workspace+"/statis_total", data)
    sys.exit(0)

##########INITIALIZATION ##########################
Factory=model.BareFactory(Map, para["Model"])
G0,W0=Factory.Build(para["Model"]["Name"], para["Lattice"]["Name"])
#IO.SaveDict("Coordinates","w", Factory.ToDict())
def Measure(G0, W0, G, W, Sigma, Polar):
    log.info("Measuring...")
    calc.Check_Denorminator(W0, Polar, Map)
    ChiTensor = calc.Calculate_ChiTensor(W0, Polar, Map)
    Chi, _ = calc.Calculate_Chi(ChiTensor, Map)
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
    data["Sigma"]=Sigma.ToDict()
    data["Polar"]=Polar.ToDict()
    data["Chi"]=Chi.ToDict()
    IO.SaveBigDict(WeightFile, data)

    #plot what you are interested in
    try:
        plot.PlotSpatial(Chi, Lat, 0, 0)
    except:
        pass

    return Chi

MessageFile=job["MessageFile"]
OutputFile=job["OutputFile"]
hist={}
hist["UnifChi"]=[]
hist["StagChi"]=[]

if len(Map.L)==2 :
    stag = Map.CoordiIndex((Map.L[0]/2, Map.L[1]/2))
elif len(Map.L)==3 :
    stag = Map.CoordiIndex((Map.L[0]/2, Map.L[1]/2, Map.L[2]/2))

if job["StartFromBare"] is True or os.path.exists(WeightFile+".pkl") is False:
    #start from bare
    Version=0
    log.info("Start from G0 and W0 to do dyson...")
    G=G0.Copy()
    W=weight.Weight("SmoothT", Map, "FourSpins", "Symmetric")
    for i in range(10):
        log.info("Round #{0}...".format(i))
        Sigma=calc.Sigma_FirstOrder(G, W, Map)
        Sigma0=calc.Sigma0_FirstOrder(G, W0, Map)
        Polar=calc.Polar_FirstOrder(G, Map)
        #######DYSON FOR W AND G###########################
        G = calc.G_Dyson(G0, Sigma0, Sigma, Map)
        W = calc.W_Dyson(W0, Polar, Map)
        ###################################################
    Chi=Measure(G0, W0, G, W, Sigma, Polar)

    Chi.FFT(1, "Space", "Time")
    hist["UnifChi"].append(Chi.Data[0,0,0,0,0,0]/Map.MaxTauBin*Map.Beta)
    hist["StagChi"].append(Chi.Data[0,0,0,0,stag,0]/Map.MaxTauBin*Map.Beta)
    IO.SaveDict(OutputFile, "w", hist)
    Chi.FFT(-1, "Space", "Time")

    parameter.BroadcastMessage(MessageFile, {"Version": Version, "Beta": Map.Beta})
    log.info("#{0} is done!".format(Version))

else:
    #########READ G,SIGMA,POLAR; CALCULATE SIGMA0 #################
    Version=parameter.GetVersion(MessageFile)
    while True:
        Version+=1
        log.info("Start #{0}...".format(Version))
        try:
            log.info("Load Sigma/Polar and G to do dyson...")
            G0,W0=Factory.Build(para["Model"]["Name"], para["Lattice"]["Name"])
            #reinitialize G0, W0 to kill accumulated error
            data=IO.LoadBigDict(WeightFile)
            G=weight.Weight("SmoothT", Map, "TwoSpins", "AntiSymmetric").FromDict(data["G"])
            paraDyson=para["Dyson"]
            MaxOrder=paraDyson["Order"]
            log.info("Collecting Sigma/Polar statistics...")
            SigmaMC, PolarMC=collect.CollectStatis(Map, MaxOrder)
            Sigma,Polar=collect.UpdateWeight(SigmaMC, PolarMC, 
                    paraDyson["ErrorThreshold"], paraDyson["OrderAccepted"])
            Sigma0=calc.Sigma0_FirstOrder(G, W0, Map)
            log.info("Dyson GW...")
            #######DYSON FOR W AND G###########################
            G = calc.G_Dyson(G0, Sigma0, Sigma, Map)
            W = calc.W_Dyson(W0, Polar,Map)
            ####### Measure ############
            Chi=Measure(G0, W0, G, W, Sigma, Polar)

            Chi.FFT(1, "Space", "Time")
            hist["UnifChi"].append(Chi.Data[0,0,0,0,0,0]/Map.MaxTauBin*Map.Beta)
            hist["StagChi"].append(Chi.Data[0,0,0,0,stag,0]/Map.MaxTauBin*Map.Beta)
            IO.SaveDict(OutputFile, "w", hist)
            Chi.FFT(-1, "Space", "Time")

            parameter.BroadcastMessage(MessageFile, {"Version": Version, "Beta": Map.Beta})
            log.info("#{0} is done!".format(Version))
        except:
            log.info("#{0} fails due to\n {1}".format(Version, traceback.format_exc()))
        finally:
            time.sleep(para["Dyson"]["SleepTime"])
