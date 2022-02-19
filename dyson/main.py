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

def Measure(para, Observable,Factory, G0, W0, G, W, SigmaDeltaT, Sigma, Polar, Determ, ChiTensor, GGGammaG=None, WWGammaW=None):
    log.info("Measuring...")
    Map=G0.Map
    ChiTensor=calc.Add_ChiTensor_ZerothOrder(ChiTensor, G, Map)
    Chi = calc.Calculate_Chi(ChiTensor, Map)

    if para["Gamma3"]:
        BKChiTensor, _ = gamma3.FullGGGammaG(GGGammaG, W0, Map)
        BKChi = gamma3.Calculate_Chi(BKChiTensor, Map)
    else:
        BKChi=Chi

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

    print "Chi=\n", np.sum(Chi.Data[0,0,0,0,:,:], axis=0)

    #print "Polar[UP,UP]=\n", Polar.Data[spinUP,0,spinUP,0,0,:]
    #print "Polar[DOWN, DOWN]=\n", Polar.Data[spinDOWN,0,spinDOWN,0,0,:]
    #print "W0=\n", W0.Data[spinUP,0,spinUP,0,1]
    # print "G0[UP,UP]=\n", G0.Data[UP,0,UP,0,0,:]
    # beta=0.5
    # Nt=16
    # t=np.array([(i+0.5)*beta/Nt for i in range(Nt)])
    # print np.exp(np.pi/2.0/beta*t*1j)/(1+1j)
    #print "G0[DOWN,DOWN]=\n", G0.Data[DOWN,0,DOWN,0,0,:]
    # G.FFT("K","T")
    # print "G[UP,UP]=\n", G.Data[UP,0,UP,0,:,-1]
    # print "G[DOWN,DOWN]=\n", G.Data[DOWN,0,DOWN,0,0,:]
    #print "SigmaDeltaT[UP,UP]=\n", SigmaDeltaT.Data[UP,0,UP,0,0]
    #print "SigmaDeltaT[DOWN,DOWN]=\n", SigmaDeltaT.Data[DOWN,0,DOWN,0,0]
    # for i in range(Map.MaxTauBin):
        # n=Sigma.Data[UP,0,UP,0,0,i]
        # print '%05f %05f %05f' % (i*Map.Beta/Map.MaxTauBin, n.real, n.imag)

    # if WWGammaW is not None:
        # RestoredWWGammaW=gamma3.UnCompressGammaW(WWGammaW, Map)
        # # RestoredWWGammaW=WWGammaW
        # print RestoredWWGammaW[0,1,0,0,:]
        # print RestoredWWGammaW[1,1,0,0,:]
        # print RestoredWWGammaW[0,0,0,0,:]
        # print RestoredWWGammaW[0,0,1,0,:]
        # print WWGammaW.shape

    data={}
    data["Chi"]=Chi.ToDict()
    data["G"]=G.ToDict()
    data["W"]=W.ToDict()
    data["W"].update(W0.ToDict())
    data["SigmaDeltaT"]=SigmaDeltaT.ToDict()
    data["Sigma"]=Sigma.ToDict()
    data["Polar"]=Polar.ToDict()
    if Map.Dim==2:
        TauSqueeze, TauRestore, TauSymFactor, TauFlag, RSqueeze, RRestore, RSymFactor=gamma3.SymmetryMapping(Map)
    if para["Gamma3"]:
        data["WWGammaW"]={"TauSqueeze": TauSqueeze.flatten().tolist(), "TauSymFactor": TauSymFactor.flatten().tolist(),
                "TauFlag": TauFlag.flatten().tolist(),
                "RSqueeze": RSqueeze.flatten().tolist(), "RSymFactor": RSymFactor.flatten().tolist(),
                "RSize": len(RRestore), "TauSize": len(TauRestore)}
        BKChiTensor.FFT("R","T")
        BKChi.FFT("R","T")
        print "BKChi=\n", np.sum(BKChi.Data[0,0,0,0,:,:], axis=0)

        data["BKChi"]=BKChi.ToDict()
        data["GGGammaG"]={"SmoothT": GGGammaG}
        if WWGammaW is not None:
            data["WWGammaW"]["SmoothT"]=gamma3.CompressGammaW(WWGammaW, Map)

    Observable.Measure(Chi, BKChi, Determ, G, Factory.NearestNeighbor)

    with DelayedInterrupt():
        try:
            log.info("Save weights into {0} File".format(WeightFile))

            IO.SaveBigDict(WeightFile, data)

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
            # plot.PlotWeightvsR("\chi", Chi,Lat,0,0)
            plot.PlotBand(G0, Lat)
        except:
            log.warning(blue("Output fails due to\n {0}".format(traceback.format_exc())))

def Dyson(IsDysonOnly, IsNewCalculation, EnforceSumRule, para, Map, Lat):
    ParaDyson=para["Dyson"]
    if not para.has_key("Version"):
        para["Version"]=0
    ########## Calulation INITIALIZATION ##########################
    Factory=model.BareFactory(Map, Lat,  para["Model"], ParaDyson["Annealing"])
    G0,W0=Factory.Build()
    # G0.FFT("K","W")
    # print G0.Data[UP,0,UP,0,0,:]
    # print "Comparison 1"
    # for n in range(Map.MaxTauBin):
        # wn=1j*(2*n+1)*np.pi/Map.Beta
        # print 1.0/(wn-4.0), G0.Data[UP,0,UP,0,0,n]*Map.Beta/Map.MaxTauBin

    # print "Comparison 2"
    # G0.FFT("K","T")
    # for t in range(Map.MaxTauBin):
        # tau=Map.IndexToTau(t)
        # G0w=-np.exp(4*tau)*(1.0-1.0/(1.0+np.exp(-Map.Beta*4)))
        # print G0.Data[UP,0,UP,0,0,t], G0w

    # print "Comparison 3"
    # for n in range(Map.MaxTauBin):
        # Gw=0.0
        # wn=(2*n+1)*np.pi/Map.Beta
        # for t in range(Map.MaxTauBin):
            # tau=Map.IndexToTau(t)
            # G0w=-np.exp(4*tau)*(1.0-1.0/(1.0+np.exp(-Map.Beta*4)))
            # Gw+=G0w*np.exp(-1j*wn*tau)*Map.Beta/Map.MaxTauBin
            # # Gw+=G0.Data[UP,0,UP,0,0,t]*np.exp(-1j*wn*tau)*Map.Beta/Map.MaxTauBin
        # print 1.0/(1j*wn-4.0), Gw

    G0.FFT("K","T")
    print G0.Data[UP,0,UP,0,0,:]
    plot.PlotTimeForList("G0UPUP_r", G0, UP, 0, UP, 0, range(Map.L[0]))
    plot.PlotBand(G0, Lat)

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
        if para["Gamma3"]:
            GGGammaG=gamma3.SimpleGG(G, Map)
    else:
        #load WeightFile, load G,W
        log.info("Load G, W from {0}".format(WeightFile))
        data=IO.LoadBigDict(WeightFile)
        G=weight.Weight("SmoothT", Map, "TwoSpins", "AntiSymmetric", "R", "T").FromDict(data["G"])
        W.FromDict(data["W"])
        SigmaDeltaT.FromDict(data["SigmaDeltaT"])
        Sigma.FromDict(data["Sigma"])
        Polar.FromDict(data["Polar"])

        if para["Gamma3"]:
            if data.has_key("GGGammaG"):
                GGGammaG=data["GGGammaG"]["SmoothT"]
                print "Read existing GGGammaG"
            else:
                GGGammaG=gamma3.SimpleGG(G, Map)

    Gold, Wold = G, W

    #while para["Version"]==0:
    while True:
        para["Version"]+=1
        log.info(green("Start Version {0}...".format(para["Version"])))
        try:
            # ratio=None   #set this will not use accumulation!
            ratio = para["Version"]/(para["Version"]+10.0)
            G0,W0=Factory.Build()
            # print W0.Data[:,0,:,0,1]
            log.info("calculating SigmaDeltaT..")
            SigmaDeltaT.Merge(ratio, calc.SigmaDeltaT_FirstOrder(G, W0, Map))
            log.info("SigmaDeltaT is done")

            # print "Polar[UP,UP]=\n", Polar.Data[spinUP,0,spinUP,0,0,:]
            # print "GammaG[UP,UP]=\n", GammaG[UP,0,:,-1]

            if IsDysonOnly or IsNewCalculation:
                log.info("accumulating Sigma/Polar statistics...")
                G = calc.G_Dyson(G0, SigmaDeltaT, Sigma, Map)
                Sigma.Merge(ratio, calc.SigmaSmoothT_FirstOrder(G, W, Map))
                log.info("calculating G...")

                G = calc.G_Dyson(G0, SigmaDeltaT, Sigma, Map)
                Polar.Merge(ratio, calc.Polar_FirstOrder(G, Map))

                if para["Gamma3"]:
                    # irreducible GGGammaG = simpleGG + GGGammaG_2 + GGGammaG_3

                    # the second term GammaG: the term from dSigma/dG
                    # GGGammaG_2 = G*(W*GGGammaG)*G
                    print "Attach W to GGGammaG"
                    GammaG = gamma3.AddW_To_GGGammaG(GGGammaG, W, G.Map)
                    print "Calculate GammaG contribution to GGGammaG"
                    GGGammaG_2 = gamma3.AddTwoG_To_GammaG(GammaG, G, G.Map)

                    if Map.MaxTauBin==Map.MaxTauBinTiny:
                        # the third term: the term from dSigma/dW
                        # GGGammaG_3 = G*((W*(G*GGGammaG)*W)*G)*G
                        print "Calculate GammaW"
                        GammaW = gamma3.AddG_To_GGGammaG(GGGammaG, G, G.Map)
                        print "Calculate WWGammaW"
                        WWGammaW = gamma3.AddTwoW_To_GammaW(GammaW, W0, W, G.Map)
                        print "Calculate GammaG from WWGammaW"
                        GammaG_FromWWGammaW=gamma3.AddG_To_WWGammaW(WWGammaW, G, G.Map)
                        print "Calculate WWGammaW contribution to GGGammaG"
                        GGGammaG_3 = gamma3.AddTwoG_To_GammaG(GammaG_FromWWGammaW, G, G.Map)
                    else:
                        WWGammaW=None
                        GGGammaG_3=0.0

                    SimpleGGGammaG=gamma3.SimpleGG(G, Map)

                    GGGammaG = SimpleGGGammaG
                    GGGammaG += +GGGammaG_2 - GGGammaG_3

            else:
                log.info("Collecting Sigma/Polar statistics...")
                SigmaStatis, PolarStatis, GammaG_MC, GammaW_MC =collect.CollectStatis(Map, para["Gamma3"])
                Sigma, Polar_MC, ParaDyson["OrderAccepted"]=collect.UpdateWeight([SigmaStatis, PolarStatis],
                        ParaDyson["ErrorThreshold"], ParaDyson["OrderAccepted"])
                #print Sigma.Data[0,0,0,0,0,0], Sigma.Data[0,0,0,0,0,-1]
                log.info("calculating G...")

                G = calc.G_Dyson(G0, SigmaDeltaT, Sigma, Map)
                SigmaDyson = calc.SigmaSmoothT_FirstOrder(G, W, Map)
                print "SigmaFromDyson=\n", SigmaDyson.Data[UP,0,UP,0,0,:]

                Polar.Merge(ratio, Polar_MC)

                if para["Gamma3"]:
                    print "Calculate WWGammaW contribution to GGGammaG"
                    WWGammaW = gamma3.AddTwoW_To_GammaW(GammaW_MC, W0, W, G.Map)
                    
                    print "Calculate GGGammaW contribution to GGGammaG"
                    GGGammaG_MC = gamma3.AddTwoG_To_GammaG(GammaG_MC, G, G.Map)
                    print "Add Simple GG contribution to GGGammaG"
                    GGGammaG = gamma3.SimpleGG(G, Map)+ GGGammaG_MC
                    # print "GammaG, mc=\n",  0.5*(np.sum(GGGammaG_MC[DOWN, :, :, :]-GGGammaG_MC[UP, :, :, :], axis=0)).diagonal()

            #######DYSON FOR W AND G###########################
            log.info("calculating W...")

            Wtmp, ChiTensor, Determ = calc.W_Dyson(W0, Polar, Map, Lat)

            if EnforceSumRule:
                ChiTensor=calc.Add_ChiTensor_ZerothOrder(ChiTensor, G, Map)
                Chi = calc.Calculate_Chi(ChiTensor, Map)
                Chi.FFT("R","T")

                while abs(Chi.Data[0,0,0,0,0,0]-0.75)>1.e-3:
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
            if para["Gamma3"]:
                Measure(para, Observable, Factory, G0, W0, G, W, SigmaDeltaT, Sigma, Polar, Determ, ChiTensor, GGGammaG, WWGammaW)
            else:
                Measure(para, Observable, Factory, G0, W0, G, W, SigmaDeltaT, Sigma, Polar, Determ, ChiTensor)
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

    para["Model"]["Symmetry"].append(para["Lattice"]["Name"])
    WeightPara={"NSublat": para["Lattice"]["NSublat"], "L":para["Lattice"]["L"],
            "Beta": float(para["Tau"]["Beta"]), "MaxTauBin": para["Tau"]["MaxTauBin"],
            "MaxTauBinTiny": para["Tau"]["MaxTauBinTiny"], "BasisNum": para["Tau"]["BasisNum"], "Symmetry": para["Model"]["Symmetry"]}
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

