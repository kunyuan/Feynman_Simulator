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
        if self.__Lat.Name in ["Square", "Cubic"]:
            Chi.FFT("K", "W")
            StagKIndex=self.__Map.CoordiIndex([e/2 for e in self.__Map.L])
            self.Append("UnifChi", Chi.Data[0,0,0,0,0,0]*Factor)
            self.Append("StagChi", Chi.Data[0,0,0,0,StagKIndex,0]*Factor)
            Chi.FFT("R","T")
            energy=np.zeros(self.__Map.MaxTauBin)
            for i in range(self.__Lat.NSublat):
                for j in range(self.__Lat.NSublat):
                    for l in NN[i][j]:
                        energy+=Chi.Data[0,i,0,j,self.__Map.CoordiIndex(l),:]/self.__Lat.NSublat
            self.Append("Energy", np.mean(energy))
        elif self.__Lat.Name in ["Checkerboard", "3DCheckerboard"]:
            Chi.FFT("K", "W")
            StagKIndex=self.__Map.CoordiIndex([0 for e in self.__Map.L])
            self.Append("UnifChi", (Chi.Data[0,0,0,0,0,0]+Chi.Data[0,0,0,1,0,0])*Factor)
            self.Append("StagChi", (Chi.Data[0,0,0,0,StagKIndex,0]-Chi.Data[0,0,0,1,StagKIndex,0])*Factor)

            Chi.FFT("R","T")
            energy=np.zeros(self.__Map.MaxTauBin)
            for i in range(self.__Lat.NSublat):
                for j in range(self.__Lat.NSublat):
                    for l in NN[i][j]:
                        energy+=Chi.Data[0,i,0,j,self.__Map.CoordiIndex(l),:]/self.__Lat.NSublat
            self.Append("Energy", np.mean(energy))

            G.FFT("R","T")
            self.Append("<S^z_A>", 0.5*(G.Data[UP,0,UP,0,0,-1]-G.Data[DOWN,0,DOWN,0,0,-1]))
            self.Append("<S^z_B>", 0.5*(G.Data[UP,1,UP,1,0,-1]-G.Data[DOWN,1,DOWN,1,0,-1]))
        elif self.__Lat.Name in ["Pyrochlore"]:
            Chi.FFT("R", "W")
            K=(4*PI, 2*PI ,0) #High symmetry point with strongest order
            self.Append("Chi_X(4Pi,2Pi,0)", 
                    self.__Lat.FourierTransformation(Chi.Data[0,:,0,:,:,0]*Factor, [K,],"Real")[1][0])
            self.Append("UnifChi", 
                    self.__Lat.FourierTransformation(Chi.Data[0,:,0,:,:,0]*Factor, [(0.0,0.0,0.0),],"Real")[1][0])
            Chi.FFT("R","T")
            energy=np.zeros(self.__Map.MaxTauBin)+1j*0
            for i in range(self.__Lat.NSublat):
                for j in range(self.__Lat.NSublat):
                    for l in NN[i][j]:
                        energy+=Chi.Data[0,i,0,j,self.__Map.CoordiIndex(l),:]/self.__Lat.NSublat
            self.Append("Energy", np.mean(energy))
        else:
            Assert(False, "model not implemented!")

        infostr="Latest measurement:\n"
        for key in self.__History.keys():
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
