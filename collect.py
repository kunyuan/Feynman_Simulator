#!/usr/bin/env python
import numpy as np
import IO 
import os, sys
sys.path.append("./calculator/") 
import lattice, logger, model, weight, parameter, plot, argparse

StatisFilePattern="_statistics"
workspace = "./"
FileList = [f for f in os.listdir(workspace) if os.path.isfile(os.path.join(workspace,f))]
StatisFileList=[f for f in FileList if f.find(StatisFilePattern) is not -1]


data={}
data["Sigma"]={}
data["Sigma"]["Histogram"]={}
data["Polar"]={}
data["Polar"]["Histogram"]={}

for quantity in {"Sigma", "Polar"}:
    i=0
    for file in StatisFileList:
        Histogram=IO.LoadBigDict(file[:-4])[quantity]["Histogram"]

        if i==0:
            data[quantity]["Histogram"] = Histogram
        else :
            data[quantity]["Histogram"]["NormAccu"] += Histogram["NormAccu"]
            data[quantity]["Histogram"]["WeightAccu"] += Histogram["WeightAccu"]
            for key in Histogram["Estimator"].keys():
                data[quantity]["Histogram"]["Estimator"][key]["Accu"] += Histogram["Estimator"][key]["Accu"] 
                data[quantity]["Histogram"]["Estimator"][key]["History"] += Histogram["Estimator"][key]["History"] 
        i+=1


IO.SaveBigDict("statistics_tot", data)

if __name__=="__main__":
    print StatisFileList

