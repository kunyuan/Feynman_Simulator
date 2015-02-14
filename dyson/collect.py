#!/usr/bin/env python
import numpy as np
import os, sys, weight
from logger import *
import parameter as para
from scipy.interpolate import LSQUnivariateSpline
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
    #define segments to do spline smoothing
    t=np.linspace(0,len(x),3)[1:-1]
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

    def MergeFromDict(self, WeightDict):
        datamat=WeightDict[self.__Weight.Name]
        """add data from another WeightEstimator"""
        self.WeightAccu+=datamat['WeightAccu']
        self.NormAccu+=datamat['NormAccu']
        if self.Norm is None:
            self.Norm=datamat['Norm']
        else:
            Assert(self.Norm==datamat['Norm'], "Norm have to be the same to merge statistics")

    def UpdateWeight(self, Name, ErrorThreshold, OrderAccepted):
        """ Weight accumulation data will be destroyed here to save memory!!!"""
        if abs(self.NormAccu)<1e-3:
            raise CollectStatisFailure("{0} 's NormAccu is 0.0!".format(Name))
        self.WeightAccu*=1.0/self.NormAccu*self.Norm
        self.NormAccu=0.0 #destroy accumulation data
        self.OrderWeight=self.WeightAccu
        MaxTauBin=self.__Map.MaxTauBin
        x=range(0, MaxTauBin)
        Shape=self.OrderWeight.shape
        Average0 = np.average(abs(np.sum(self.OrderWeight[:,0,0,0,0,0,:], axis=0)))
        Info=[]
        DimList=[]
        for order in xrange(1,Shape[0]+1):
            # Order Index is equal to real Order-1
            RelativeError=0.0
            weight=self.OrderWeight[order-1,...]
            for index, _ in np.ndenumerate(weight[...,0]):
                sp1, sub1, sp2, sub2, vol=index
                y=weight[index] # y is a function of tau
                if not np.allclose(y, 0.0, 1e-5): 
                    smooth, sigma=Smooth(x, y)
                    relative=abs(sigma)/abs(Average0)
                    if relative>RelativeError:
                        RelativeError=relative
                        error=sigma
                        Original=weight[index].copy()
                        Smoothed=smooth.copy()
                        Position=(order,sp1,sub1,sp2,sub2,vol)
                    #weight[index]=smooth
            log.info("Maximum at Order {0} is {1}".format(order, RelativeError))
            IsAccpted=RelativeError<ErrorThreshold or order<=OrderAccepted
            State="Accepted with relative error {0:.2g}".format(RelativeError, ErrorThreshold)
            if not IsAccpted:
                State="NOT "+State
            Info.append([Name, x, Original, Smoothed, error, Position, State])
            if not IsAccpted:
                order-=1
                break
        try:
            self.__Plot(Info)
        except:
            raise
            log.warning("Failed to plot statistics of {0} at order {1}".format(Name, Position[0]))
        NewOrderAccepted=order
        log.info("OrderAccepted={0}".format(NewOrderAccepted))

        Dict={self.__Weight.Name: np.sum(self.OrderWeight[:NewOrderAccepted], axis=0)}
        self.__Weight.FromDict(Dict)
        return self.__Weight, NewOrderAccepted

    def __Plot(self, Info):
        import matplotlib.pyplot as plt
        color=plt.cm.rainbow(np.linspace(0,1,len(Info)))
        fig=plt.figure()
        ax1=plt.subplot(1, 2, 1)
        ax2=plt.subplot(1, 2, 2)
        InfoStr=""
        for i in range(len(Info)):
            Name, x, y, smooth, sigma, Position, State=Info[i]
            mid=len(x)/2
            ax1.plot(x, y.real, '+', c=color[i])
            ax1.plot(x, smooth.real, '-k')
            ax1.errorbar(x[mid], smooth[mid].real, yerr=sigma.real,
                    label="Error {0:.2g}".format(sigma.real), c=color[i],elinewidth=3)
            ax1.set_xlim([x[0],x[-1]])
            ax2.plot(x, y.imag, '+', c=color[i],
                    label="Order {0}/Sp:{1},Sub:{2},Coord:{3}".format(Position[0], 
                     (Position[1], Position[3]), (Position[2],Position[4]),
                     self.__Map.IndexToCoordi(Position[5])))
            ax2.plot(x, smooth.imag, '-k')
            ax2.errorbar(x[mid], smooth[mid].imag, yerr=sigma.imag, 
                    label="Error {0:.2g}\n{1}".format(sigma.imag, State), c=color[i],elinewidth=3)
            ax2.set_xlim([x[0],x[-1]])

        ax1.legend(loc='best', fancybox=True, framealpha=0.5, prop={'size':6})
        ax2.legend(loc='best', fancybox=True, framealpha=0.5, prop={'size':6})
        ax1.set_xlabel("$\\tau_{bin}$")
        ax2.set_xlabel("$\\tau_{bin}$")
        plt.savefig("{0}_Smoothed.pdf".format(Name))
        plt.close()

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
            SigmaSmoothT.MergeFromDict(Dict['Sigma']['Histogram'])
            PolarSmoothT.MergeFromDict(Dict['Polar']['Histogram'])
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
    Sigma,SigmaOrder=SigmaSmoothT.UpdateWeight("Sigma", ErrorThreshold, OrderAccepted["Sigma"])
    Polar,PolarOrder=PolarSmoothT.UpdateWeight("Polar", ErrorThreshold, OrderAccepted["Polar"])
    log.info("Accepted Sigma order : {0}; Accepted Polar order : {1}".
            format(SigmaOrder, PolarOrder))
    if SigmaOrder==0 or PolarOrder==0:
        raise CollectStatisFailure("Either Sigma or Polar's OrderAccepted is still zero!")
    OrderAccepted={"Sigma": SigmaOrder, "Polar": PolarOrder}
    return Sigma, Polar, OrderAccepted

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
