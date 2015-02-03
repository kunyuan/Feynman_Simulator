#!/usr/bin/env python
import numpy as np
import calculator as calc
import lattice as lat
import collect
from weight import UP,DOWN,IN,OUT
from logger import *
import os, sys, model, weight, measure, parameter, plot, argparse, time, traceback, signal
import plot

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--PID", help="use PID to find the input file")
parser.add_argument("-f", "--file", help="use file path to find the input file")
parser.add_argument("-c", "--collect", action="store_true", help="collect all the _statis.* file into statis_total.*")

########## Environment INITIALIZATION ##########################
args = parser.parse_args()
if args.PID:
    InputFile=os.path.join(workspace, "infile/_in_DYSON_"+str(args.PID))
elif args.file:
    InputFile=os.path.abspath(args.file)
else:
    Assert(False, "Do not understand the argument!")

job, para=parameter.Load(InputFile)
ParaFile="{0}_DYSON_para".format(job["PID"])
WeightFile=job["WeightFile"]
MessageFile=job["MessageFile"]
OutputFile=job["OutputFile"]
StatisFile=os.path.join(workspace, "statis_total")

DoesWeightFileExist=os.path.exists(WeightFile+".hkl")
if DoesWeightFileExist:
    try:
        log.info("Load previous DYSHON_para file")
        para=IO.LoadDict(ParaFile)["Para"]
    except:
        log.info("Previous DYSHON_para file, use _in_DYSON_ file as para instead")

WeightPara={"NSublat": para["Lattice"]["NSublat"], "L":para["Lattice"]["L"],
            "Beta": float(para["Tau"]["Beta"]), "MaxTauBin": para["Tau"]["MaxTauBin"]}
ParaDyson=para["Dyson"]
Map=weight.IndexMap(**WeightPara)
Lat=lat.Lattice(para["Lattice"]["Name"], Map)

if args.collect:
    log.info("Collect statistics only...")
    SigmaMC, PolarMC=collect.CollectStatis(Map, ParaDysob["Order"])
    data ={}
    data["Sigma"] = {"Histogram": SigmaMC.ToDict()}
    data["Polar"] = {"Histogram": PolarMC.ToDict()}
    with DelayedInterrupt():
        IO.SaveBigDict(StatisFile, data)
    sys.exit(0)

########## Calulation INITIALIZATION ##########################
Factory=model.BareFactory(Map, Lat,  para["Model"], para["Dyson"]["Annealing"])
G0,W0=Factory.Build()
IO.SaveDict("Coordinates","w", Factory.ToDict())

Observable=measure.Observable(Map, Lat)

def Measure(G0, W0, G, W, SigmaDeltaT, Sigma, Polar, Determ, ChiTensor):
    log.info("Measuring...")
    Chi = calc.Calculate_Chi(ChiTensor, Map)

    ##########OUTPUT AND FILE SAVE ####################
    spinUP=Map.Spin2Index(UP,UP)
    spinDOWN=Map.Spin2Index(DOWN,DOWN)

    Polar.FFT("R","T")
    W0.FFT("R")
    W.FFT("R","T")
    G.FFT("R","T")
    SigmaDeltaT.FFT("R")
    Sigma.FFT("R","T")
    Chi.FFT("R","T")

    #print "Polar=\n", Polar.Data[spinUP,0,spinUP,0,0,:]
    #print "W=\n", W.Data[spinUP,0,spinUP,0,0,:]
    #print "G[UP,UP]=\n", G.Data[UP,0,UP,0,0,:]
    #print "G[DOWN,DOWN]=\n", G.Data[UP,0,UP,0,0,:]
    #print "SigmaDeltaT=\n", SigmaDeltaT.Data[UP,0,UP,0,0]
    #print "Sigma=\n", Sigma.Data[UP,0,UP,0,0,:]
    #print "Chi=\n", Chi.Data[0,0,0,0,1,:]

    data={}
    data["Chi"]=Chi.ToDict()
    data["G"]=G.ToDict()
    data["W"]=W.ToDict()
    data["W"].update(W0.ToDict())
    data["SigmaDeltaT"]=SigmaDeltaT.ToDict()
    data["Sigma"]=Sigma.ToDict()
    data["Polar"]=Polar.ToDict()
    Observable.Measure(Chi, Determ, G)

    with DelayedInterrupt():
        log.info("Save weights into {0} File".format(WeightFile))
        IO.SaveBigDict(WeightFile, data)
        parameter.Save(ParaFile, para)  #Save Parameters
        Observable.Save(OutputFile)

    #plot what you are interested in
    try:
        plot.PlotChi(Chi, Lat)
    except:
        pass

#if MessageFile does not exist, Version will be 0
Version=parameter.GetVersion(MessageFile)

if DoesWeightFileExist is True:
    #load WeightFile, load G,W
    log.info("Load G, W from {0}".format(WeightFile))
    data=IO.LoadBigDict(WeightFile)
    G=weight.Weight("SmoothT", Map, "TwoSpins", "AntiSymmetric", "R", "T").FromDict(data["G"])
    W=weight.Weight("SmoothT", Map, "FourSpins", "Symmetric","R","T").FromDict(data["W"])
else:
    #not load WeightFile
    log.info("Start from bare G, W")
    G=G0.Copy()
    W=weight.Weight("SmoothT", Map, "FourSpins", "Symmetric","R","T")

Gold, Wold = G, W
SigmaDeltaT=weight.Weight("DeltaT", Map, "TwoSpins", "AntiSymmetric","R")

if (job["DysonOnly"] is True) or (DoesWeightFileExist is False):
    #not load StatisFile
    Sigma=weight.Weight("SmoothT", Map, "TwoSpins", "AntiSymmetric","R","T")
    Polar=weight.Weight("SmoothT", Map, "FourSpins", "Symmetric","R","T")

#while True:
while Version<30:
    Version+=1
    log.info("Start Version {0}...".format(Version))
    try:
        ratio = Version/(Version+10.0)
        G0,W0=Factory.Build()
        SigmaDeltaT.Merge(ratio, calc.SigmaDeltaT_FirstOrder(G, W0, Map))

        if (job["DysonOnly"] is True) or (DoesWeightFileExist is False):
            log.info("accumulating Sigma/Polar statistics...")
            Sigma.Merge(ratio, calc.SigmaSmoothT_FirstOrder(G, W, Map))
            Polar.Merge(ratio, calc.Polar_FirstOrder(G, Map))
        else:
            log.info("Collecting Sigma/Polar statistics...")
            MaxOrder=ParaDyson["Order"]
            SigmaMC, PolarMC=collect.CollectStatis(Map, MaxOrder)
            Sigma, Polar, SigmaOrder, PolarOrder=collect.UpdateWeight(SigmaMC, PolarMC, 
                    ParaDyson["ErrorThreshold"], ParaDyson["OrderAccepted"])
            if SigmaOrder==0 or PolarOrder==0:
                log.info("Version {0} fails due to Sigma or Polar not accepted. \n".format(Version))
                continue
        try:
            #######DYSON FOR W AND G###########################
            print "calculating W..."
            W, ChiTensor, Determ = calc.W_Dyson(W0, Polar, Map)
            print "calculating G..."
            G = calc.G_Dyson(G0, SigmaDeltaT, Sigma, Map)
        except KeyboardInterrupt, SystemExit:
            raise
        except:
            Factory.RevertField(ParaDyson["Annealing"])
            G, W = Gold, Wold
        else:
            Gold, Wold = G, W
            Measure(G0, W0, G, W, SigmaDeltaT, Sigma, Polar, Determ, ChiTensor)
            Factory.DecreaseField(ParaDyson["Annealing"])

        log.info("Version {0} is done!".format(Version))
        parameter.BroadcastMessage(MessageFile, {"Version": Version, "Beta": Map.Beta})
    except KeyboardInterrupt, SystemExit:
        log.info("Terminating Dyson\n {1}".format(Version, traceback.format_exc()))
        sys.exit(0)
    except:
        log.info("Version {0} fails due to\n {1}".format(Version, traceback.format_exc()))
    finally:
        if (job["DysonOnly"] is False) and (DoesWeightFileExist is True):
            time.sleep(ParaDyson["SleepTime"])
