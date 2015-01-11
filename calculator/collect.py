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
        if Name=="SmoothT":
            self.Beta=self.Map.Beta
            self.Shape.append(self.Map.MaxTauBin)
            self.__HasTau=True
        elif Name=="DeltaT":
            self.__HasTau=False
        else:
            Assert(False, "Should be either .SmoothT or .DeltaT in Name, not {0}".format(Name))
        self.Data=np.zeros(self.Shape, dtype=complex)
    def FromDict(self, data):
        log.info("Loading {0} Matrix...".format(self.Name));
        if self.Name in data:
            log.info("Load {0}".format(self.Name))
            datamat = data[self.Name]
            ######RESHAPE data[self.Name]
            OldShape=[self.NSpin**2, self.NSublat**2]+self.__OriginShape[VOL:]
            MidShape=[self.NSpin, self.NSpin, self.NSublat,self.NSublat]+self.__OriginShape[VOL:]
            NewShape=self.__OriginShape
            self.__AssertShape(datamat.shape, OldShape)
            datamat=datamat.reshape(MidShape).swapaxes(1,2).reshape(NewShape)
            self.__AssertShape(self.Shape, datamat.shape)
            self.Data=datamat
        else:
            Assert(False, "{0} not found!").format(self.Name)


def GetFileList():
    FileList = [f for f in os.listdir(workspace) if os.path.isfile(os.path.join(workspace,f))]
    StatisFileList=[f for f in FileList if f.find(StatisFilePattern) is not -1]
    return StatisFileList

def CollectStatis(_map, _order):
    _FileList=GetFileList()
    log.info("Collect statistics from {0}".format(_FileList))
    SigmaSmoothT=WeightEstimator("SmoothT", _map, "TwoSpins", _order)
    PolarSmoothT=WeightEstimator("SmoothT", _map, "FourSpins", _order)

if __name__=="__main__":
    WeightPara={"NSublat": 1, "L":[16, 16],
                "Beta": 0.1, "MaxTauBin": 64}
    map=weight.IndexMap(**WeightPara)
    Order=3
    CollectStatis(map, Order)
