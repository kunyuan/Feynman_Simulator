#!/usr/bin/env python
import numpy as np
import os, sys, weight
from logger import *
import parameter as para

StatisFilePattern="_statis"

class WeightEstimator():
    def __init__(self, Weight, Order): 
        self.Shape=[Order, Weight.NSpin**2, Weight.NSublat**2]+Weight.Shape[Weight.VOLDIM:]
        self.__Weight=Weight
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
            print self.Norm, _WeightEstimator.Norm
            Assert(self.Norm==_WeightEstimator.Norm, "Norm have to be the same to merge statistics")
    
    def GetWeight(self, ErrorThreshold, OrderAccepted):
        """add all different orders together"""
        Dict={self.__Weight.Name: np.sum(self.WeightAccu, axis=0)/
                           self.NormAccu*self.Norm*self.Map.MaxTauBin/self.Map.Beta}
        self.__Weight.FromDict(Dict)
        return self.__Weight

    def FromDict(self, data):
        datamat=data[self.__Weight.Name]
        self.WeightAccu=datamat['WeightAccu']
        self.NormAccu=datamat['NormAccu']
        self.Norm=datamat['Norm']
        self.__AssertShape(self.Shape, self.WeightAccu.shape)

    def ToDict(self):
        datatmp={}
        datatmp["WeightAccu"]=self.WeightAccu
        datatmp["NormAccu"]=self.NormAccu
        datatmp["Norm"]=self.Norm
        return {self.__Weight.Name: datatmp}

    def __AssertShape(self, shape1, shape2):
        Assert(tuple(shape1)==tuple(shape2), \
                "Shape {0} is expected instead of shape {1}!".format(shape1, shape2))

def GetFileList():
    FileList = [f for f in os.listdir(workspace) if os.path.isfile(os.path.join(workspace,f))]
    StatisFileList=[os.path.join(workspace, f) for f in FileList if f.find(StatisFilePattern) is not -1]
    return StatisFileList

def CollectStatis(_map, _order):
    Sigma=weight.Weight("SmoothT", _map, "TwoSpins", "AntiSymmetric")
    SigmaSmoothT=WeightEstimator(Sigma, _order)
    Polar=weight.Weight("SmoothT", _map, "FourSpins", "Symmetric")
    PolarSmoothT=WeightEstimator(Polar, _order)
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

def UpdateWeight(SigmaSmoothT, PolarSmoothT, ErrorThreshold, OrderAccepted):
    Sigma=SigmaSmoothT.GetWeight(ErrorThreshold, OrderAccepted)
    Polar=PolarSmoothT.GetWeight(ErrorThreshold, OrderAccepted)
    return Sigma, Polar

if __name__=="__main__":
    WeightPara={"NSublat": 1, "L":[4, 4],
                "Beta": 0.5, "MaxTauBin":128}
    Order=3
    map=weight.IndexMap(**WeightPara)

    SigmaSmoothT, PolarSmoothT=CollectStatis(map, Order)

    data ={}
    data["Sigma"] = {"Histogram": SigmaSmoothT.ToDict()}
    data["Polar"] = {"Histogram": PolarSmoothT.ToDict()}

    IO.SaveBigDict(workspace+"/statis_total", data)
