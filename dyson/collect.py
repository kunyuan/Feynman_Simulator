#!/usr/bin/env python
import numpy as np
import os, sys, weight
from logger import *
import parameter as para
import traceback

StatisFilePattern="_statis"
AcceptRatio=0.75

class CollectStatisFailure(Exception):
    def __init__(self, msg):
       self.Message = msg
    def __str__(self):
        return str(self.Message)

def Smooth(x,y):
    """return: smoothed funciton s(x) and estimation of sigma of y for one data point"""
    from scipy.interpolate import LSQUnivariateSpline
    #define segments to do spline smoothing
    t=np.linspace(0,len(x),10)[1:-1]
    sr = LSQUnivariateSpline(x, y.real, t)
    si = LSQUnivariateSpline(x, y.imag, t)
    return sr(x)+si(x)*1j, (sr.get_residual()/len(x))**0.5+(si.get_residual()/len(x))**0.5*1j

class WeightEstimator():
    def __init__(self, Weight, Order): 
        self.Shape=[Order,]+Weight.Shape
        self.__Map=Weight.Map
        self.__Weight=Weight
        self.WeightAccu=np.zeros(self.Shape, dtype=complex)
        self.NormAccu=0.0
        self.Norm=None
        self.Order=Order

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
    
    def GetNewOrderAccepted(self, Name, ErrorThreshold, OrderAccepted):
        self.OrderWeight=self.WeightAccu/self.NormAccu*self.Norm
        if abs(self.NormAccu)<1e-3:
            raise CollectStatisFailure("{0} 's NormAccu is 0.0!".format(Name))
        MaxTauBin=self.__Map.MaxTauBin
        x=range(0, MaxTauBin)
        Shape=self.OrderWeight.shape
        Average0 = np.average(abs(np.sum(self.OrderWeight[:,0,0,0,0,0,:], axis=0)))
        Info=[]
        for orderindex in range(0, Shape[0]):
            # Order Index is equal to real Order-1
            RelativeError=0
            for sp1 in range(Shape[1]):
                for sub1 in range(Shape[2]):
                    for sp2 in range(Shape[3]):
                        for sub2 in range(Shape[4]):
                            for vol in range(Shape[5]):
                                y=self.OrderWeight[orderindex,sp1,sub1,sp2,sub2,vol, :]
                                if not np.allclose(y, 0.0, 1e-10): 
                                    smooth, sigma=Smooth(x, y)
                                    relative=abs(sigma)/abs(Average0)
                                    if abs(relative)>RelativeError:
                                        RelativeError=relative
                                        error=sigma
                                        Original=y.copy()
                                        Smoothed=smooth.copy()
                                        Position=(orderindex,(sp1,sp2),(sub1,sub2),vol)
                                    self.OrderWeight[orderindex,sp1,sub1,sp2,sub2,vol, :]=smooth
            State="Accepted with relative error {0:.2g}".format(RelativeError, ErrorThreshold)
            if RelativeError>=ErrorThreshold:
                State="NOT "+State
            Info.append([Name, x, Original, Smoothed, error, Position, State])
            log.info("Maximum at Order {0} is {1}".format(orderindex, RelativeError))
            if RelativeError>=ErrorThreshold:
                orderindex -= 1
                break
        try:
            self.__Plot(Info)
        except:
            log.warning("Failed to plot statistics of {0} at order {1}".format(Name, Position[0]))
        if orderindex+1>OrderAccepted:
            NewOrderAccepted=orderindex+1
        else:
            NewOrderAccepted=OrderAccepted
        log.info("OrderAccepted={0}".format(NewOrderAccepted))
        return NewOrderAccepted

    def GetWeight(self, OrderAccepted):
        OrderIndex=OrderAccepted-1
        Dict={self.__Weight.Name: np.sum(self.OrderWeight[:OrderIndex+1], axis=0)}
        self.__Weight.FromDict(Dict)
        return self.__Weight

    def __Plot(self, Info):
        import matplotlib.pyplot as plt
        fig=plt.figure()
        ax1=plt.subplot(1, 2, 1)
        ax2=plt.subplot(1, 2, 2)
        InfoStr=""
        for Name, x, y, smooth, sigma, Position, State in Info:
            mid=len(x)/2
            Order=Position[0]+1
            ax1.plot(x, y.real, '+', x, smooth.real, '-')
            ax1.errorbar(x[mid], smooth[mid].real, yerr=sigma.real,
                    label="Error {0:.2g}".format(sigma.real))
            ax1.set_xlim([x[0],x[-1]])
            ax2.plot(x, y.imag, '+', label="Order {0}/Sp:{1},Sub:{2},Coord:{3}".format(Order, 
                     Position[1], Position[2], self.__Map.IndexToCoordi(Position[3])))
            ax2.plot(x, smooth.imag, '-')
            ax2.errorbar(x[mid], smooth[mid].imag, yerr=sigma.imag, 
                    label="Error {0:.2g}\n{1}".format(sigma.imag, State))
            ax2.set_xlim([x[0],x[-1]])

        ax1.legend(loc='best', fancybox=True, framealpha=0.5, prop={'size':6})
        ax2.legend(loc='best', fancybox=True, framealpha=0.5, prop={'size':6})
        ax1.set_xlabel("$\\tau_{bin}$")
        ax2.set_xlabel("$\\tau_{bin}$")
        plt.savefig("{0}_Smoothed.pdf".format(Name))
        plt.close()

    def FromDict(self, data):
        datamat=data[self.__Weight.Name]
        self.WeightAccu=datamat['WeightAccu']
        self.NormAccu=datamat['NormAccu']
        self.Norm=datamat['Norm']
        self.OrderWeight=self.WeightAccu*self.NormAccu/self.Norm
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
    FileList = [f for f in FileList if f[0]!="_"]
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
    if len(_FileList)==0:
        raise CollectStatisFailure("No statistics files to read!") 
    log.info("Collect statistics from {0}".format(_FileList))
    Total=len(_FileList)
    Success=0.0
    for f in _FileList:
        try:
            log.info("Merging {0} ...".format(f));
            Dict=IO.LoadBigDict(f)
            SigmaTemp.FromDict(Dict['Sigma']['Histogram'])
            PolarTemp.FromDict(Dict['Polar']['Histogram'])
            SigmaSmoothT.Merge(SigmaTemp)
            PolarSmoothT.Merge(PolarTemp)
        except:
            log.info("Fails to merge\n {0}".format(traceback.format_exc()))
        else:
            Success+=1.0
    log.info("{0}/{1} statistics files read!".format(int(Success), Total))
    if float(Success)/Total<AcceptRatio:
        raise CollectStatisFailure("More than {0}% statistics files fail to read!".format(100.0*AcceptRatio)) 
    return (SigmaSmoothT, PolarSmoothT)

def UpdateWeight(StatisCollected, ErrorThreshold, OrderAccepted):
    SigmaSmoothT, PolarSmoothT=StatisCollected
    SigmaOrder=SigmaSmoothT.GetNewOrderAccepted("Sigma", ErrorThreshold, OrderAccepted)
    PolarOrder=PolarSmoothT.GetNewOrderAccepted("Polar", ErrorThreshold, OrderAccepted)
    log.info("Accepted Sigma order : {0}; Accepted Polar order : {1}".
            format(SigmaOrder, PolarOrder))
    Sigma=SigmaSmoothT.GetWeight(SigmaOrder)
    Polar=PolarSmoothT.GetWeight(PolarOrder)
    if SigmaOrder==0 or PolarOrder==0:
        raise CollectStatisFailure("Either Sigma or Polar's OrderAccepted is still zero!")
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
