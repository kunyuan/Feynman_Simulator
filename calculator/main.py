#!/usr/bin/env python
import numpy as np
import calculator as calc
import lattice as lat
import collect
from weight import UP,DOWN,IN,OUT
from logger import *
import os, sys, model, weight, measure, parameter, plot, argparse, time, traceback
import plot, gc

def debug_log_common_objects():
    print "Logging common objects"
    import gc
    print "garbage collector got", gc.collect()
    import objgraph
 
    # Reimplement show_most_common_types, but write to file
    import operator
    stats = sorted(objgraph.typestats().items(), key=operator.itemgetter(1),
                   reverse=True)
    limit = 30
    if limit:
        stats = stats[:limit]
    width = max(len(name) for name, count in stats)
    out=["#START %s" % time.time()] + ["%s %s" % (name.ljust(width), count) for name, count in stats[:limit]]
    f=open("most_common_object.txt", "a")
    f.write("\n".join(out) + "\n")
    f.close()
    print "Logging common objects - done"

def start_pdb(signal, trace):
    import pdb
    pdb.set_trace()

import signal
signal.signal(signal.SIGINT, start_pdb)

#gc.enable()
#gc.set_debug(gc.DEBUG_LEAK)
# show the dirt ;-)
#dump_garbage()
#gc.collect()

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
    InputFile=os.path.join(workspace, "infile/_in_DYSON_0")

    #Assert(False, "Do not understand the argument!")

job, para=parameter.Load(InputFile)
ParaFile="{0}_DYSON_para".format(job["PID"])
WeightFile=job["WeightFile"]
MessageFile=job["MessageFile"]
OutputFile=job["OutputFile"]
StatisFile=os.path.join(workspace, "statis_total")

DoesWeightFileExist=os.path.exists(WeightFile+".hkl")
if DoesWeightFileExist:
    try:
        log.info(green("Load previous DYSHON_para file"))
        para=parameter.LoadPara(ParaFile)
    except:
        log.warning(red("Previous DYSHON_para file, use _in_DYSON_ file as para instead"))

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

def Measure(Version, G0, W0, G, W, SigmaDeltaT, Sigma, Polar, Determ, ChiTensor):
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

    #print "Polar[UP,UP]=\n", Polar.Data[spinUP,0,spinUP,0,0,:]
    #print "Polar[DOWN, DOWN]=\n", Polar.Data[spinDOWN,0,spinDOWN,0,0,:]
    #print "W=\n", W.Data[spinUP,0,spinUP,0,0,:]
    #print "G[UP,UP]=\n", G.Data[UP,0,UP,0,0,:]
    #print "G[DOWN,DOWN]=\n", G.Data[DOWN,0,DOWN,0,0,:]
    #print "SigmaDeltaT[UP,UP]=\n", SigmaDeltaT.Data[UP,0,UP,0,0]
    #print "SigmaDeltaT[DOWN,DOWN]=\n", SigmaDeltaT.Data[DOWN,0,DOWN,0,0]
    #print "Sigma=\n", Sigma.Data[UP,0,UP,0,0,:]
    #print "Chi=\n", Chi.Data[0,0,0,0,1,:]
    #print "Chi=\n", Chi.Data[0,0,0,0,:,0]

    data={}
    data["Chi"]=Chi.ToDict()
    data["G"]=G.ToDict()
    data["W"]=W.ToDict()
    data["W"].update(W0.ToDict())
    data["SigmaDeltaT"]=SigmaDeltaT.ToDict()
    data["Sigma"]=Sigma.ToDict()
    data["Polar"]=Polar.ToDict()
    #Observable.Measure(Chi, Determ, G)

    with DelayedInterrupt():
        try:
            log.info("Save weights into {0} File".format(WeightFile))
            IO.SaveBigDict(WeightFile, data)
            para["Version"]=Version
            parameter.Save(ParaFile, para)  #Save Parameters
            Observable.Save(OutputFile)
            #plot what you are interested in
            #plot.PlotChi(Chi, Lat)
            #plot.PlotSpatial(Chi, Lat, 0, 0, 0) 
            #plot.PlotChi_2D(Chi, Lat)
        except:
            log.info("Output fails due to\n {0}".format(traceback.format_exc()))

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

if not para.has_key("Version"):
    Version=0
else:
    Version=para["Version"]

while True:
#while Version<1:
    Version+=1
    log.info(green("Start Version {0}...".format(Version)))
    try:
        ratio = Version/(Version+10.0)
        G0,W0=Factory.Build()
        SigmaDeltaT.Merge(ratio, calc.SigmaDeltaT_FirstOrder(G, W0, Map))
        #SigmaDeltaT=calc.SigmaDeltaT_FirstOrder(G, W0, Map)

        if job["DysonOnly"] or (DoesWeightFileExist is False):
            log.info("accumulating Sigma/Polar statistics...")
            #Sigma.Merge(ratio, calc.SigmaSmoothT_FirstOrder(G, W, Map))
            Sigma=calc.SigmaSmoothT_FirstOrder(G, W, Map)
            log.info("calculating G...")
            G = calc.G_Dyson(G0, SigmaDeltaT, Sigma, Map)
            #Polar.Merge(ratio, calc.Polar_FirstOrder(G, Map))
            Polar=calc.Polar_FirstOrder(G, Map)
        else:
            log.info("Collecting Sigma/Polar statistics...")
            Statis=collect.CollectStatis(Map, ParaDyson["Order"])
            Sigma, Polar=collect.UpdateWeight(Statis, ParaDyson["ErrorThreshold"], ParaDyson["OrderAccepted"])
            log.info("calculating G...")
            G = calc.G_Dyson(G0, SigmaDeltaT, Sigma, Map)
        #######DYSON FOR W AND G###########################
        log.info("calculating W...")
        W, ChiTensor, Determ = calc.W_Dyson(W0, Polar, Map, Lat)

    except calc.DenorminatorTouchZero as err:
        #failure due to denorminator touch zero
        Factory.RevertField(ParaDyson["Annealing"])
        log.warning(green("Version {0} fails due to:\n{1}".format(Version,err)))
        G, W = Gold, Wold
    except collect.CollectStatisFailure as err:
        #failure due to statis files collection
        log.warning(green("Version {0} fails due to:\n{1}".format(Version,err)))
        G, W = Gold, Wold
    except KeyboardInterrupt, SystemExit:
        #exit
        log.info("Terminating Dyson\n {1}".format(Version, traceback.format_exc()))
        sys.exit(0)
    except:
        #unknown reason failure, just fail dyson completely for safty
        log.error(red("Dyson fails due to\n {1}".format(Version, traceback.format_exc())))
        G, W = Gold, Wold
        sys.exit(0)
    else:
        #everything works prefectly 
        Gold, Wold = G, W
        Measure(Version, G0, W0, G, W, SigmaDeltaT, Sigma, Polar, Determ, ChiTensor)
        Factory.DecreaseField(ParaDyson["Annealing"])
        log.info("Version {0} is done!".format(Version))
        parameter.BroadcastMessage(MessageFile, {"Version": Version, "Beta": Map.Beta})
    finally:
        if (job["DysonOnly"] is False) and (DoesWeightFileExist is True):
            time.sleep(ParaDyson["SleepTime"])
        gc.collect()
        gc.collect()
        debug_log_common_objects()



