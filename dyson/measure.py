import numpy as np
import calculator as calc
import lattice as lat
from weight import UP,DOWN,IN,OUT
from logger import *
import os, sys, weight
PI=np.pi

class Observable:
    def __init__(self, Map, lat):
        self.__History={}
        self.__Lat=lat
        self.__Map=Map

    def Append(self, Name, Value):
        if self.__History.has_key(Name) is False:
            self.__History[Name]=[]
        self.__History[Name].append(Value)

    def Measure(self, Chi, Determinate, G, NN):
        self.Append("1-JP", Determinate.min())
        Factor=self.__Map.Beta/self.__Map.MaxTauBin
        if self.__Lat.Name in ["Square", "Cubic", "Checkerboard", "3DCheckerboard"]:
            Chi.FFT("R", "W")
            UnifK=(0.0,)*self.__Map.Dim
            StagK=(PI,)*self.__Map.Dim
            _, ChiK=self.__Lat.FourierTransformation(Chi.Data[0,:,0,:,:,0]*Factor, [UnifK,StagK],"Real")
            self.Append("UnifChi", ChiK[0])
            self.Append("StagChi", ChiK[1])
            self.Append("StagMag^2 density", ChiK[1]/self.__Map.Vol/self.__Map.NSublat/self.__Map.Beta)

        elif self.__Lat.Name in ["Pyrochlore"]:
            Chi.FFT("R", "W")
            KList=[(0.0,0.0,0.0), (4*PI, 2*PI ,0)] #High symmetry point with strongest order
            _, ChiK=self.__Lat.FourierTransformation(Chi.Data[0,:,0,:,:,0]*Factor, KList,"Real")
            self.Append("UnifChi", ChiK[0])
            self.Append("Chi_X(4Pi,2Pi,0)", ChiK[1])
        else:
            Assert(False, "model not implemented!")

        Chi.FFT("R","W")
        energy=0j
        for i in range(self.__Lat.NSublat):
            for j in range(self.__Lat.NSublat):
                for l in NN[i][j]:
                    energy+=Chi.Data[0,i,0,j,self.__Map.CoordiIndex(l),0]/self.__Lat.NSublat
        self.Append("Energy", energy/self.__Map.MaxTauBin)

        G.FFT("R","T")
        for i in range(self.__Map.NSublat):
            self.Append("<Sz_{0}>".format(i), 0.5*(G.Data[UP,i,UP,i,0,-1]-G.Data[DOWN,i,DOWN,i,0,-1]))

        infostr="Latest measurement:\n"
        for key in sorted(self.__History.keys()):
            infostr+="{0}={1}\n".format(key, self.__History[key][-1])
        log.info(infostr)

    def Load(self, FileName):
        try:
            Dict=IO.LoadDict(FileName)
        except:
            Dict=self.__History  #use old History if we fail to load
        finally:
            self.__History=Dict

    def Save(self, FileName):
        #TODO: add error bar estimation
        IO.SaveDict(FileName, "w", self.__History)
