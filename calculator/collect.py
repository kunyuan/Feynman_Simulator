#!/usr/bin/env python
import numpy as np
import os, sys, weight
from logger import *
import parameter as para

StatisFilePattern="_statis"

class WeightEstimator():
    def __init__(self, Name, Map, NSpin, Order): 
        """Name: 'SmoothT' or 'DeltaT'
           NSpin: 'TwoSpins' or 'FourSpins'
        """
        self.Map=Map
        if NSpin is "TwoSpins":
            self.NSpin=2
        elif NSpin is "FourSpins":
            self.NSpin=4
        else:
            Assert(False, "Only accept TwoSpins or FourSpins, not {0}".format(NSpin))
        self.Shape=[Order, self.NSpin**2, self.Map.NSublat**2, self.Map.Vol]
        self.Name=Name
        if Name=="SmoothT":
            self.Beta=self.Map.Beta
            self.Shape.append(self.Map.MaxTauBin)
            self.__HasTau=True
        elif Name=="DeltaT":
            self.__HasTau=False
        else:
            Assert(False, "Should be either .SmoothT or .DeltaT in Name, not {0}".format(Name))
        self.WeightAccu=np.zeros(self.Shape, dtype=complex)
        self.NormAccu=0.0
        self.Norm=None

    def Copy(self):
        """return a deep copy of Weight instance"""
        import copy
        return copy.deepcopy(self)

    def Merge(self, _WeightEstimator):
        """add data from another WeightEstimator"""
        self.WeightAccu+=_WeightEstimator.WeightAccu
        self.NormAccu+=_WeightEstimator.NormAccu
        if self.Norm is None:
            self.Norm=_WeightEstimator.Norm
        else:
            Assert(self.Norm==_WeightEstimator.Norm, "Norm have to be the same to merge statistics")
    
    def GetWeight(self, ErrorThreshold, OrderAccepted):
        """add all different orders together"""
        return np.sum(self.WeightAccu, axis=0)/self.NormAccu*self.Norm*self.Map.MaxTauBin/self.Map.Beta

    def FromDict(self, data):
        datamat=data[self.Name]
        self.WeightAccu=datamat['WeightAccu']
        self.NormAccu=datamat['NormAccu']
        self.Norm=datamat['Norm']
        self.__AssertShape(self.Shape, self.WeightAccu.shape)

    def ToDict(self):
        datatmp={}
        datatmp["WeightAccu"]=self.WeightAccu
        datatmp["NormAccu"]=self.NormAccu
        datatmp["Norm"]=self.Norm
        return {self.Name: datatmp}

    def __AssertShape(self, shape1, shape2):
        Assert(tuple(shape1)==tuple(shape2), \
                "Shape {0} is expected instead of shape {1}!".format(shape1, shape2))

def GetFileList():
    FileList = [f for f in os.listdir(workspace) if os.path.isfile(os.path.join(workspace,f))]
    StatisFileList=[os.path.join(workspace, f) for f in FileList if f.find(StatisFilePattern) is not -1]
    return StatisFileList

def CollectStatis(_map, _order):
    SigmaSmoothT=WeightEstimator("SmoothT", _map, "TwoSpins", _order)
    PolarSmoothT=WeightEstimator("SmoothT", _map, "FourSpins", _order)
    SigmaTemp=SigmaSmoothT.Copy()
    PolarTemp=PolarSmoothT.Copy()
    _FileList=GetFileList()
    log.info("Collect statistics from {0}".format(_FileList))
    for f in _FileList:
        log.info("Merging {0} ...".format(f));
        Dict=IO.LoadBigDict(f)
        SigmaTemp.FromDict(Dict['Sigma']['Histogram'])
        PolarTemp.FromDict(Dict['Polar']['Histogram'])
        SigmaSmoothT.Merge(SigmaTemp)
        PolarSmoothT.Merge(PolarTemp)
    return (SigmaSmoothT, PolarSmoothT)

def UpdateWeight(_map, _order, ErrorThreshold, OrderAccepted, SigmaSmoothT, PolarSmoothT):
    Sigma=weight.Weight("SmoothT", _map, "TwoSpins", "AntiSymmetric")
    Sigma.Data=SigmaSmoothT.GetWeight(ErrorThreshold, OrderAccepted)
    Polar=weight.Weight("SmoothT", _map, "FourSpins", "Symmetric")
    Polar.Data=PolarSmoothT.GetWeight(ErrorThreshold, OrderAccepted)

if __name__=="__main__":
    WeightPara={"NSublat": 1, "L":[4, 4],
                "Beta": 0.5, "MaxTauBin":128}
    Order=4
    map=weight.IndexMap(**WeightPara)

    SigmaSmoothT, PolarSmoothT=CollectStatis(map, Order)

    data ={}
    data["Sigma"] = {"Histogram": SigmaSmoothT.ToDict()}
    data["Polar"] = {"Histogram": PolarSmoothT.ToDict()}

    IO.SaveBigDict(workspace+"/statis_total", data)
