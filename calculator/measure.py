import numpy as np
import calculator as calc
import lattice as lat
from weight import UP,DOWN,IN,OUT
from logger import *
import os, sys, weight

class Observable:
    def __init__(self, Map, lat):
        self.__History={}
        self.__Lat=lat
        self.__Map=Map

    def Append(self, Name, Value):
        if self.__History.has_key(Name) is False:
            self.__History[Name]=[]
        self.__History[Name].append(Value)

    def Measure(self, Chi, Determinate, G):
        self.Append("1-JP", Determinate.min())
        Chi.FFT("K", "W")
        Factor=self.__Map.Beta/self.__Map.MaxTauBin
        if self.__Lat.Name in ["Square", "Cubic"]:
            StagKIndex=self.__Map.CoordiIndex([e/2 for e in self.__Map.L])
            self.Append("UnifChi", Chi.Data[0,0,0,0,0,0]*Factor)
            self.Append("StagChi", Chi.Data[0,0,0,0,StagKIndex,0]*Factor)
        elif self.__Lat.Name in ["Checkerboard", "3DCheckerboard"]:
            StagKIndex=self.__Map.CoordiIndex([0 for e in self.__Map.L])
            self.Append("UnifChi", (Chi.Data[0,0,0,0,0,0]+Chi.Data[0,0,0,1,0,0])*Factor)
            self.Append("StagChi", (Chi.Data[0,0,0,0,StagKIndex,0]-Chi.Data[0,0,0,1,StagKIndex,0])*Factor)
            G.FFT("R","T")
            self.Append("<S^z_A>", 0.5*(G.Data[UP,0,UP,0,0,-1]-G.Data[DOWN,0,DOWN,0,0,-1]))
            self.Append("<S^z_B>", 0.5*(G.Data[UP,1,UP,1,0,-1]-G.Data[DOWN,1,DOWN,1,0,-1]))
        elif self.__Lat.Name in ["Pyrochlore"]:
            return
            #StagKIndex=self.__Map.CoordiIndex([0 for e in self.__Map.L])
            #self.__History["UnifChi"].append((Chi.Data[0,0,0,0,0,0]+Chi.Data[0,0,0,1,0,0])*Factor)
            #self.__History["StagChi"].append((Chi.Data[0,0,0,0,StagKIndex,0]-Chi.Data[0,0,0,1,StagKIndex,0])*Factor)
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
