#!/usr/bin/env python
import numpy as np
import calculator as calc
import lattice as lat
import collect
from weight import UP,DOWN,IN,OUT
from logger import *
import os, sys, model, weight, measure, parameter, plot, argparse, time, traceback, signal

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

#def handler(signal, frame):
    #print "ctrl-c received"
    #sys.exit(0)
#signal.signal(signal.SIGINT, handler)

job, para=parameter.Load(InputFile)
ParaFile="{0}_DYSON_para".format(job["PID"])
WeightFile=job["WeightFile"]
MessageFile=job["MessageFile"]
OutputFile=job["OutputFile"]
StatisFile=os.path.join(workspace, "statis_total")

DoesWeightFileExist=os.path.exists(WeightFile+".hkl")
if DoesWeightFileExist:
    try:
        para=IO.LoadDict(ParaFile)["Para"]
    except:
        pass

WeightPara={"NSublat": para["Lattice"]["NSublat"], "L":para["Lattice"]["L"],
            "Beta": para["Tau"]["Beta"], "MaxTauBin": para["Tau"]["MaxTauBin"]}
ParaDyson=para["Dyson"]
Map=weight.IndexMap(**WeightPara)
Lat=lat.Lattice(para["Lattice"]["Name"], Map)

if args.collect:
    MaxOrder=para["Dyson"]["Order"]
    SigmaMC, PolarMC=collect.CollectStatis(Map, MaxOrder)
    data ={}
    data["Sigma"] = {"Histogram": SigmaMC.ToDict()}
    data["Polar"] = {"Histogram": PolarMC.ToDict()}
    IO.SaveBigDict(StatisFile, data)
    sys.exit(0)

########## Calulation INITIALIZATION ##########################
Factory=model.BareFactory(Map, Lat,  para["Model"], para["Dyson"]["Annealing"])
G0,W0=Factory.Build()
IO.SaveDict("Coordinates","w", Factory.ToDict())

Observable=measure.Observable(Map, Lat)

def Measure(G0, W0, G, W, Sigma0, Sigma, Polar, Determ, ChiTensor):
    log.info("Measuring...")
    Chi = calc.Calculate_Chi(ChiTensor, Map)

    ##########OUTPUT AND FILE SAVE ####################
    spinUP=Map.Spin2Index(UP,UP)
    spinDOWN=Map.Spin2Index(DOWN,DOWN)

    Polar.FFT("R","T")
    W0.FFT("R")
    W.FFT("R","T")
    G.FFT("R","T")
    Sigma0.FFT("R")
    Sigma.FFT("R","T")
    Chi.FFT("R","T")

    #print "Polar=\n", Polar.Data[spinUP,0,spinUP,0,0,:]
    #print "W=\n", W.Data[spinUP,0,spinUP,0,0,:]
    #print "G[UP,UP]=\n", G.Data[UP,0,UP,0,0,:]
    #print "G[DOWN,DOWN]=\n", G.Data[UP,0,UP,0,0,:]
    #print "Sigma0=\n", Sigma0.Data[UP,0,UP,0,0]
    #print "Sigma=\n", Sigma.Data[UP,0,UP,0,0,:]
    #print "Chi=\n", Chi.Data[0,0,0,0,1,:]

    data={}
    data["Chi"]=Chi.ToDict()
    data["G"]=G.ToDict()
    data["W"]=W.ToDict()
    data["W"].update(W0.ToDict())
    data["Sigma0"]=Sigma0.ToDict()
    data["Sigma"]=Sigma.ToDict()
    data["Polar"]=Polar.ToDict()
    Observable.Measure(Chi, Determ)

    with DelayedKeyboardInterrupt():
        log.info("Save weights into {0} File".format(WeightFile))
        IO.SaveBigDict(WeightFile, data)
        parameter.Save(ParaFile, para)  #Save Parameters
        Observable.Save(OutputFile)

    #plot what you are interested in

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
Sigma0=weight.Weight("DeltaT", Map, "TwoSpins", "AntiSymmetric","R","T")

if job["DysonOnly"] is True:
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
        Sigma0.Merge(ratio, calc.SigmaDeltaT_FirstOrder(G, W0, Map))

        if job["DysonOnly"] is True:
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
            G = calc.G_Dyson(G0, Sigma0, Sigma, Map)
        except KeyboardInterrupt, SystemExit:
            raise
        except:
            Factory.RevertField(ParaDyson["Annealing"])
            G, W = Gold, Wold
        else:
            Gold, Wold = G, W
            Measure(G0, W0, G, W, Sigma0, Sigma, Polar, Determ, ChiTensor)
            Factory.DecreaseField(ParaDyson["Annealing"])
        log.info("Version {0} is done!".format(Version))
        parameter.BroadcastMessage(MessageFile, {"Version": Version, "Beta": Map.Beta})
    except KeyboardInterrupt, SystemExit:
        log.info("Terminating Dyson\n {1}".format(Version, traceback.format_exc()))
        sys.exit(0)
    except:
        log.info("Version {0} fails due to\n {1}".format(Version, traceback.format_exc()))
    finally:
        if job["DysonOnly"] is False:
            time.sleep(ParaDyson["SleepTime"])

#if job["StartFromBare"] is True or os.path.exists(WeightFile+".hkl") is False:
    ##start from bare
    #Version=0
    #log.info("Start from G0 and W0 to do dyson...")
    #G=G0.Copy()
    #W=weight.Weight("SmoothT", Map, "FourSpins", "Symmetric","R","T")
    #Gold, Wold = G, W
    #Sigma=weight.Weight("SmoothT", Map, "TwoSpins", "AntiSymmetric","R","T")
    #Sigma0=weight.Weight("DeltaT", Map, "TwoSpins", "AntiSymmetric","R","T")
    #Polar=weight.Weight("SmoothT", Map, "FourSpins", "Symmetric","R","T")
    #for i in range(30):
        #ratio = i/(i+10.0)
        #log.info("Round #{0}...".format(i))
        #G0,W0=Factory.Build(para["Model"]["Name"], para["Lattice"]["Name"])

        #Sigma.Merge(ratio, calc.SigmaSmoothT_FirstOrder(G, W, Map))
        #Sigma0.Merge(ratio, calc.SigmaDeltaT_FirstOrder(G, W0, Map))
        #Polar.Merge(ratio, calc.Polar_FirstOrder(G, Map))

        #try:
            ########DYSON FOR W AND G###########################
            #print "calculating W..."
            #W, ChiTensor, Determ = calc.W_Dyson(W0, Polar, Map)
            #print "calculating G..."
            #G = calc.G_Dyson(G0, Sigma0, Sigma, Map)
        #except:
            #Factory.RevertField(para["Dyson"]["Annealing"])
            #G, W = Gold, Wold
        #else:
            #Gold, Wold = G, W
            #Measure(G0, W0, G, W, Sigma0, Sigma, Polar, Determ, ChiTensor)
            #Factory.DecreaseField(para["Dyson"]["Annealing"])

        ####################################################

    #parameter.BroadcastMessage(MessageFile, {"Version": Version, "Beta": Map.Beta})
    #log.info("Version {0} is done!".format(Version))

#else:
    ##########READ G,SIGMA,POLAR; CALCULATE SIGMA0 #################
    #Version=parameter.GetVersion(MessageFile)
    #Observable.Load(OutputFile)
    #data=IO.LoadBigDict(WeightFile)
    #G=weight.Weight("SmoothT", Map, "TwoSpins", "AntiSymmetric").FromDict(data["G"])

    #while True:
        #Version+=1
        #log.info("Start #{0}...".format(Version))
        #try:
            #log.info("Load Sigma/Polar and G to do dyson...")
            #G0,W0=Factory.Build(para["Model"]["Name"], para["Lattice"]["Name"])
            ##reinitialize G0, W0 to kill accumulated error, and change with externalfield
            #paraDyson=para["Dyson"]
            #MaxOrder=paraDyson["Order"]
            #log.info("Collecting Sigma/Polar statistics...")
            #SigmaMC, PolarMC=collect.CollectStatis(Map, MaxOrder)

            #Sigma, Polar, SigmaOrder, PolarOrder=collect.UpdateWeight(SigmaMC, PolarMC, 
                    #paraDyson["ErrorThreshold"], paraDyson["OrderAccepted"])

            #if SigmaOrder==0 or PolarOrder==0:
                #log.info("#{0} fails due to Sigma or Polar not accepted. \n".format(Version))
                #continue
                
            #Sigma0=calc.SigmaDeltaT_FirstOrder(G, W0, Map)
            #log.info("Dyson GW...")

            #try:
                ########DYSON FOR W AND G###########################
                #print "calculating W..."
                #W, ChiTensor, Determ = calc.W_Dyson(W0, Polar, Map)
                #print "calculating G..."
                #G = calc.G_Dyson(G0, Sigma0, Sigma, Map)
            #except:
                #Factory.RevertField(para["Dyson"]["Annealing"])
                #G, W = Gold, Wold
            #else:
                #Gold, Wold = G, W
                #Measure(G0, W0, G, W, Sigma0, Sigma, Polar, Determ, ChiTensor)
                #Factory.DecreaseField(para["Dyson"]["Annealing"])
            #log.info("Version {0} is done!".format(Version))
            #parameter.BroadcastMessage(MessageFile, {"Version": Version, "Beta": Map.Beta})

        #except:
            #log.info("#{0} fails due to\n {1}".format(Version, traceback.format_exc()))
        #finally:
            #time.sleep(para["Dyson"]["SleepTime"])
