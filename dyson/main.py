#!/usr/bin/env python
from logger import *
import sys
#log.info("Python Version: {0}".format(sys.version))
import numpy as np
import calculator as calc
import gamma3_calc as gamma3
import lattice as lat
import collect
from weight import UP,DOWN,IN,OUT
import os, model, weight, measure, parameter, plot, argparse, time, traceback
import plot, gc
#def start_pdb(signal, trace):
    #import pdb
    #pdb.set_trace()

#import signal
#start in pdb mode after Ctrl-C
#signal.signal(signal.SIGINT, start_pdb)

def Measure(para, Observable,Factory, G0, W0, G, W, SigmaDeltaT, Sigma, Polar, Determ, ChiTensor, GammaG, GammaW):
    log.info("Measuring...")
    ChiTensor=calc.Add_ChiTensor_ZerothOrder(ChiTensor, G, Map)
    Chi = calc.Calculate_Chi(ChiTensor, Map)

    ##########OUTPUT AND FILE SAVE ####################
    spinUP=Map.Spin2Index(UP,UP)
    spinDOWN=Map.Spin2Index(DOWN,DOWN)

    Polar.FFT("R","T")
    W0.FFT("R")
    W.FFT("R","T")
    G0.FFT("R","T")
    G.FFT("R","T")
    SigmaDeltaT.FFT("R")
    Sigma.FFT("R","T")
    Chi.FFT("R","T")

    #print "Polar[UP,UP]=\n", Polar.Data[spinUP,0,spinUP,0,0,:]
    #print "Polar[DOWN, DOWN]=\n", Polar.Data[spinDOWN,0,spinDOWN,0,0,:]
    #print "W0=\n", W0.Data[spinUP,0,spinUP,0,1]
    # print "G0[UP,UP]=\n", G0.Data[UP,0,UP,0,0,:]
    # beta=0.5
    # Nt=16
    # t=np.array([(i+0.5)*beta/Nt for i in range(Nt)])
    # print np.exp(np.pi/2.0/beta*t*1j)/(1+1j)
    #print "G0[DOWN,DOWN]=\n", G0.Data[DOWN,0,DOWN,0,0,:]
    #print "G[UP,UP]=\n", G.Data[UP,0,UP,0,0,:]
    #print "G[DOWN,DOWN]=\n", G.Data[DOWN,0,DOWN,0,0,:]
    #print "SigmaDeltaT[UP,UP]=\n", SigmaDeltaT.Data[UP,0,UP,0,0]
    #print "SigmaDeltaT[DOWN,DOWN]=\n", SigmaDeltaT.Data[DOWN,0,DOWN,0,0]
    # for i in range(Map.MaxTauBin):
        # n=Sigma.Data[UP,0,UP,0,0,i]
        # print '%05f %05f %05f' % (i*Map.Beta/Map.MaxTauBin, n.real, n.imag)

    print "Chi, mc=\n",  0.5*(np.sum(GammaG[DOWN, :, :, :]-GammaG[UP, :, :, :], axis=0)).diagonal()

    # GammaG_simple = calc.SimpleGG(G, Map)
    # GGW=calc.GGW(GammaG_simple, W, G.Map)
    # GGWGG=calc.AddTwoGToGammaG(GGW, G, G.Map)
    # print "Chi, dyson=\n",  0.5*(np.sum(GGWGG[DOWN, :, :, :]-GGWGG[UP, :, :, :], axis=0)).diagonal()

    # GGW=calc.GGW(GammaG, W, G.Map)
    # GGWGG=calc.AddTwoGToGammaG(GGW, G, G.Map)

    # GGammaG = calc.AddG_To_GammaG(GammaG, G, G.Map)
    # WWGammaW=calc.WWGammaW(GGammaG, W0, W, G.Map)
    # GammaGFromGammaW=calc.GammaWToGammaG(WWGammaW, G, G.Map)

    print "GammaW, type0, diagonal, mc=\n", GammaW[0, 1, 1, :, :].diagonal()
    print "GammaW, type4, diagonal, mc=\n", GammaW[4, 1, 1, :, :].diagonal()

    print "GammaW, type0, tau1=0,  mc=\n",  GammaW[0, 1, 1, 0, :]
    print "GammaW, type4, tau1=0,  mc=\n",  GammaW[4, 1, 1, 0, :]
    print "GammaW, type5, tau1=0,  mc=\n",  GammaW[5, 1, 1, 0, :]

    print "GammaG, UP, r=0, mc=\n", GammaG[UP, 0, :, :].diagonal()
    print "GammaG, UP, r=1, mc=\n", GammaG[UP, 1, :, :].diagonal()

    # GammaGFirstOrder=calc.GammaG_FirstOrder(GammaG, G, W0, Map)
    # SimpleGammaG=calc.SimpleGG(G, Map)

    # GammaG_dyson =SimpleGammaG+GammaGFirstOrder
    # GammaG_dyson += +GGWGG - GammaGFromGammaW

    # print "Chi, dyson=\n",  0.5*(np.sum(GammaG_dyson[DOWN, :, :, :]-GammaG_dyson[UP, :, :, :], axis=0)).diagonal()

    #_,ChiTensor,_=calc.W_Dyson(W0, Polar, Polar.Map, Lat) 
    #ChiTensor.FFT("R","T")
    #print "ChiTensor=\n", np.sum(ChiTensor.Data[spinUP,0,spinUP,0,:,:], axis=0)
    #print "ChiTensor=\n", np.sum(ChiTensor.Data[spinUP,0,spinDOWN,0,:,:], axis=0)
    print "Chi, polar=\n", np.sum(Chi.Data[0,0,0,0,:,:], axis=0)

    # print "WWGammaW, type0, diagonal, mc=\n", GammaW[0, 1, 1, :, :].diagonal()
    # print "WWGammaW, type0, diagonal, dyson=\n", WWGammaW[0, 1, 1, :, :].diagonal()

    # print "WWGammaW, type4, diagonal, mc=\n", GammaW[4, 1, 1, :, :].diagonal()
    # print "WWGammaW, type4, diagonal, dyson=\n", WWGammaW[4, 1, 1, :, :].diagonal()

    # GammaG_simple = calc.SimpleGG(G, G.Map)
    # GGammaG = calc.AddG_To_GammaG(GammaG_simple, G, G.Map)
    # print "GGammaG, type0, dyson=\n", GGammaG[0, 0, 0, :, :].diagonal()
    # WWGammaW_dyson=calc.WWGammaW(GGammaG, W0, W, G.Map)
    # print "WWGammaW, type0, dyson, diagonal=\n", WWGammaW_dyson[0, 1, 1, :, :].diagonal()
    # print "WWGammaW, type4, dyson, t1=0,=\n", WWGammaW_dyson[4, 1, 1, 0, :]

    # GammaG_dyson = calc.GammaWToGammaG(WWGammaW_dyson, G, G.Map)
    # print "GammaG, UP, dyson=\n", GammaG_dyson[UP, 1, :, :].diagonal()

    # print "GammaG, last term, DOWN, dyson=\n", -GammaG_dyson[DOWN, 1, :, :].diagonal()
    # print "GammaG, DOWN, mc=\n", GammaG[DOWN, 1, :, :].diagonal()

    # print "WWGammaW, type4, dyson=\n", np.sum(WWGammaW_dyson[4, 1, 1, :, :])/WWGammaW_dyson.shape[3]/WWGammaW_dyson.shape[4]
    # print "WWGammaW, type4, dyson, diagonal=\n", WWGammaW_dyson[4, 1, 1, 0, :]

    # print "WWGammaW, type5, dyson=\n", np.sum(WWGammaW_dyson[5, 1, 1, :, :])/WWGammaW_dyson.shape[3]/WWGammaW_dyson.shape[4]
    # print "WWGammaW, type5, dyson, diagonal=\n", WWGammaW_dyson[5, 1, 1, 0, :]

    data={}
    data["Chi"]=Chi.ToDict()
    data["G"]=G.ToDict()
    data["W"]=W.ToDict()
    data["W"].update(W0.ToDict())
    data["SigmaDeltaT"]=SigmaDeltaT.ToDict()
    data["Sigma"]=Sigma.ToDict()
    data["Polar"]=Polar.ToDict()
    data["GammaG"]={"SmoothT": GammaG}
    if GammaW is not None:
        data["GammaW"]={"SmoothT": GammaW}
    Observable.Measure(Chi, Determ, G, Factory.NearestNeighbor)

    with DelayedInterrupt():
        try:
            log.info("Save weights into {0} File".format(WeightFile))

            #####TODO: NOT UPDATING WEIGHT FILE
            IO.SaveBigDict(WeightFile, data)
            #####################################

            parameter.Save(ParaFile, para)  #Save Parameters
            Observable.Save(OutputFile)

            #plot what you are interested in
            plot.PlotChiAlongPath(Chi, Lat)
            plot.PlotTime("G", G, UP, 0, UP, 0, 0)
            plot.PlotTime("G0UPUP", G0, UP, 0, UP, 0, 0)
            plot.PlotTime("G0DOWNDOWN", G0, DOWN, 0, DOWN, 0, 0)
            plot.PlotTime("Sigma", Sigma, DOWN, 0, DOWN, 0, 0)
            plot.PlotTime("Polar", Polar, spinUP, 0, spinUP, 0, 0)
            plot.PlotSpatial(Chi, Lat, 0, 0) 
            plot.PlotChi_2D(Chi, Lat)
            plot.PlotWeightvsR("\chi", Chi,Lat,0,0)
        except:
            log.info(blue("Output fails due to\n {0}".format(traceback.format_exc())))

def Dyson(IsDysonOnly, IsNewCalculation, EnforceSumRule, para, Map, Lat):
    ParaDyson=para["Dyson"]
    if not para.has_key("Version"):
        para["Version"]=0
    ########## Calulation INITIALIZATION ##########################
    Factory=model.BareFactory(Map, Lat,  para["Model"], ParaDyson["Annealing"])
    G0,W0=Factory.Build()
    IO.SaveDict("Coordinates","w", Factory.ToDict())
    Observable=measure.Observable(Map, Lat)
    W=weight.Weight("SmoothT", Map, "FourSpins", "Symmetric","R","T")
    SigmaDeltaT=weight.Weight("DeltaT", Map, "TwoSpins", "AntiSymmetric","R")
    Sigma=weight.Weight("SmoothT", Map, "TwoSpins", "AntiSymmetric","R","T")
    Polar=weight.Weight("SmoothT", Map, "FourSpins", "Symmetric","R","T")

    if IsNewCalculation:
        #not load WeightFile
        log.info("Start from bare G, W")
        G=G0.Copy()
        GammaG=gamma3.SimpleGG(G, Map)
        GammaW=np.zeros([6, Map.Vol, Map.Vol, Map.MaxTauBin, Map.MaxTauBin])+0.0*1j
    else:
        #load WeightFile, load G,W
        log.info("Load G, W from {0}".format(WeightFile))
        data=IO.LoadBigDict(WeightFile)
        G=weight.Weight("SmoothT", Map, "TwoSpins", "AntiSymmetric", "R", "T").FromDict(data["G"])
        W.FromDict(data["W"])
        SigmaDeltaT.FromDict(data["SigmaDeltaT"])
        Sigma.FromDict(data["Sigma"])
        Polar.FromDict(data["Polar"])

        if data.has_key("GammaG"):
            GammaG=data["GammaG"]["SmoothT"]
            print "Read existing GammaG"
        else:
            GammaG=gamma3.SimpleGG(G, Map)

        if data.has_key("GammaW"):
            GammaW=data["GammaW"]["SmoothT"]
            print "Read existing GammaW"
        else:
            GammaW=np.zeros([6, Map.Vol, Map.Vol, Map.MaxTauBin, Map.MaxTauBin])+0.0*1j

    Gold, Wold = G, W

    #while para["Version"]==0:
    while True:
        para["Version"]+=1
        log.info(green("Start Version {0}...".format(para["Version"])))
        try:
            ratio=None   #set this will not use accumulation!
            #ratio = para["Version"]/(para["Version"]+10.0)
            G0,W0=Factory.Build()
            # print W0.Data[:,0,:,0,1]
            log.info("calculating SigmaDeltaT..")
            SigmaDeltaT.Merge(ratio, calc.SigmaDeltaT_FirstOrder(G, W0, Map))
            log.info("SigmaDeltaT is done")

            # print "Polar[UP,UP]=\n", Polar.Data[spinUP,0,spinUP,0,0,:]
            # print "GammaG[UP,UP]=\n", GammaG[UP,0,:,-1]

            if IsDysonOnly or IsNewCalculation:
                log.info("accumulating Sigma/Polar statistics...")
                Sigma.Merge(ratio, calc.SigmaSmoothT_FirstOrder(G, W, Map))
                log.info("calculating G...")

                G = calc.G_Dyson(G0, SigmaDeltaT, Sigma, Map)
                Polar.Merge(ratio, calc.Polar_FirstOrder(G, Map))

                GGW=gamma3.GGW(GammaG, W, G.Map)
                GGWGG=gamma3.AddTwoGToGammaG(GGW, G, G.Map)

                GGammaG = gamma3.AddG_To_GammaG(GammaG, G, G.Map)
                # WWGammaW1=gamma3.WWGammaW(GGammaG, W0, W, G.Map)
                WWGammaW=gamma3.FastWWGammaW(GGammaG, W0, W, G.Map)
                GammaGFromGammaW=gamma3.GammaWToGammaG(WWGammaW, G, G.Map)

                # GammaGFirstOrder=gamma3.GammaG_FirstOrder(GammaG, G, W0, Map)
                GammaGFirstOrder=gamma3.FastGammaG_RPA(GammaG, G, W0, Map)
                SimpleGammaG=gamma3.SimpleGG(G, Map)

                GammaG =SimpleGammaG+GammaGFirstOrder
                GammaG += +GGWGG - GammaGFromGammaW

                print "GammaGFirstOrder, dyson=\n",  0.5*(np.sum(GammaGFirstOrder[DOWN, :, :, :]-GammaGFirstOrder[UP, :, :, :], axis=0)).diagonal()
                # print "FastGammaGRPA, dyson=\n",  0.5*(np.sum(FastGammaGRPA[DOWN, :, :, :]-FastGammaGRPA[UP, :, :, :], axis=0)).diagonal()
                # print "GGWGG, mc=\n",  0.5*(np.sum(GGWGG[DOWN, :, :, :]-GGWGG[UP, :, :, :], axis=0)).diagonal()
                #print "GGWGG, mc=\n",  GGWGG[UP, 0, :, :].diagonal()
                # print "GammaGFromGammaW, mc=\n",  0.5*(np.sum(GammaGFromGammaW[DOWN, :, :, :]-GammaGFromGammaW[UP, :, :, :], axis=0)).diagonal()
                print "SimpleGammaG, dyson=\n",  0.5*(np.sum(SimpleGammaG[DOWN, :, :, :]-SimpleGammaG[UP, :, :, :], axis=0)).diagonal()
                print "GammaG, dyson=\n",  0.5*(np.sum(GammaG[DOWN, :, :, :]-GammaG[UP, :, :, :], axis=0)).diagonal()
                print "WWGammaW, type0, diagonal, dyson=\n", WWGammaW[0, 1, 1, 0, :]
                print "WWGammaW, type4, diagonal, dyson=\n", WWGammaW[4, 1, 1, 0, :]
                print "WWGammaW, type5, diagonal, dyson=\n", WWGammaW[5, 1, 1, 0, :]

            else:
                log.info("Collecting Sigma/Polar statistics...")
                SigmaStatis, PolarStatis, GammaG_MC, GammaW_MC=collect.CollectStatis(Map)
                Sigma, Polar, ParaDyson["OrderAccepted"]=collect.UpdateWeight([SigmaStatis, PolarStatis],
                        ParaDyson["ErrorThreshold"], ParaDyson["OrderAccepted"])
                #print Sigma.Data[0,0,0,0,0,0], Sigma.Data[0,0,0,0,0,-1]
                log.info("calculating G...")

                G = calc.G_Dyson(G0, SigmaDeltaT, Sigma, Map)
                SigmaDyson = calc.SigmaSmoothT_FirstOrder(G, W, Map)
                print "SigmaFromDyson=\n", SigmaDyson.Data[UP,0,UP,0,0,:]

                GammaW = gamma3.WWGammaW(GammaW_MC, W0, W, G.Map)
                GGGammaG_MC = gamma3.AddTwoGToGammaG(GammaG_MC, G, G.Map)
                GammaG = gamma3.SimpleGG(G, Map)+ gamma3.GammaG_FirstOrder(GammaG, G, W0, Map) + GGGammaG_MC
                print "GammaG, mc=\n",  0.5*(np.sum(GGGammaG_MC[DOWN, :, :, :]-GGGammaG_MC[UP, :, :, :], axis=0)).diagonal()

            #######DYSON FOR W AND G###########################
            log.info("calculating W...")
            Wtmp, ChiTensor, Determ = calc.W_Dyson(W0, Polar, Map, Lat)
            if EnforceSumRule:
                ChiTensor=calc.Add_ChiTensor_ZerothOrder(ChiTensor, G, Map)
                Chi = calc.Calculate_Chi(ChiTensor, Map)
                Chi.FFT("R","T")

                while abs(Chi.Data[0,0,0,0,0,0]-0.75)>1.e-2:
                    SumRuleRatio = np.sqrt(0.75/Chi.Data[0,0,0,0,0,0])
                    PolarSumRule = Polar
                    PolarSumRule.Data = PolarSumRule.Data*SumRuleRatio
                    Wtmp, ChiTensor, Determ = calc.W_Dyson(W0, Polar, Map, Lat)
                    ChiTensor=calc.Add_ChiTensor_ZerothOrder(ChiTensor, G, Map)
                    Chi = calc.Calculate_Chi(ChiTensor, Map)
                    Chi.FFT("R","T")
                    print "Chi(r=0,t=0)", Chi.Data[0,0,0,0,0,0]
            W = Wtmp


        except calc.DenorminatorTouchZero as err:
            #failure due to denorminator touch zero
            log.info(green("Version {0} fails due to:\n{1}".format(para["Version"],err)))
            Factory.RevertField(ParaDyson["Annealing"])
            G, W = Gold, Wold
            SigmaDeltaT.RollBack()
            Sigma.RollBack()
            Polar.RollBack()

        except collect.CollectStatisFailure as err:
            #failure due to statis files collection
            log.info(green("Version {0} fails due to:\n{1}".format(para["Version"],err)))
            G, W = Gold, Wold
            SigmaDeltaT.RollBack()
            Sigma.RollBack()
            Polar.RollBack()
        except KeyboardInterrupt, SystemExit:
            #exit
            log.info("Terminating Dyson\n {1}".format(para["Version"], traceback.format_exc()))
            sys.exit(0)
        except:
            #unknown reason failure, just fail dyson completely for safty
            log.info(red("Dyson fails due to\n {1}".format(para["Version"], traceback.format_exc())))
            sys.exit(0)
        else:
            #everything works prefectly 
	    log.info("everything is going well!")
            Gold, Wold = G, W
            Measure(para, Observable, Factory, G0, W0, G, W, SigmaDeltaT, Sigma, Polar, Determ, ChiTensor, GammaG, GammaW)
            IsSuccessed=Factory.DecreaseField(ParaDyson["Annealing"])
            # Factor=2.0 if IsSuccessed else 1.0
            # Factor=2.0
            Factor = 1.0 + 3.0/(1.0*para["Version"]+1.0)
            parameter.BroadcastMessage(MessageFile, 
                    {"Version": para["Version"], "Beta": Map.Beta, "SqueezeFactor": Factor})
            log.info("Version {0} is done!".format(para["Version"]))
        finally:
            log.info(green("Memory Usage before collecting: {0} MB".format(memory_usage())))
            gc.collect()
            log.info(green("Memory Usage : {0} MB".format(memory_usage())))
            if not IsDysonOnly and not IsNewCalculation:
                time.sleep(ParaDyson["SleepTime"])

if __name__=="__main__":
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
    global ParaFile
    ParaFile="{0}_DYSON_para".format(job["PID"])
    global WeightFile
    WeightFile=job["WeightFile"]
    global MessageFile
    MessageFile=job["MessageFile"]
    global OutputFile
    OutputFile=job["OutputFile"]
    global StatisFile
    StatisFile=os.path.join(workspace, "statis_total")

    IsNewCalculation=not os.path.exists(WeightFile+".hkl")
    if not IsNewCalculation: 
        try:
            log.info(green("Try to load previous DYSHON_para file"))
            para=parameter.LoadPara(ParaFile)
            log.info("{0} para file is loaded".format(ParaFile))
        except:
            log.warning(red("Previous DYSHON_para file does not exist, use _in_DYSON_ file as para instead"))
            parameter.Save(ParaFile, para)  #Save Parameters

    WeightPara={"NSublat": para["Lattice"]["NSublat"], "L":para["Lattice"]["L"],
                "Beta": float(para["Tau"]["Beta"]), "MaxTauBin": para["Tau"]["MaxTauBin"]}
    Map=weight.IndexMap(**WeightPara)
    Lat=lat.Lattice(para["Lattice"]["Name"], Map)

    if args.collect:
        log.info("Collect statistics only...")
        SigmaMC, PolarMC=collect.CollectStatis(Map)
        collect.UpdateWeight((SigmaMC, PolarMC), para["Dyson"]["ErrorThreshold"], para["Dyson"]["OrderAccepted"])
        data ={}
        data["Sigma"] = {"Histogram": SigmaMC.ToDict()}
        data["Polar"] = {"Histogram": PolarMC.ToDict()}
        with DelayedInterrupt():
            IO.SaveBigDict(StatisFile, data)
        sys.exit(0)
    else:
        Dyson(job["DysonOnly"],  IsNewCalculation, job["SumRule"], para, Map, Lat)

    log.info("calculation ended!")

