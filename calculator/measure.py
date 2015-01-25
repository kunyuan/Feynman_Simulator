import numpy as np
import calculator as calc
import lattice as lat
from logger import *
import os, sys, weight

class Observable:
    def __init__(self, Map, lat):
        self.__History={"UnifChi":[],"StagChi":[],"1-JP":[]}
        self.__Lat=lat
        self.__Map=Map

    def Measure(self, Chi, Determinate):
        self.__History["1-JP"].append(Determinate.min())
        Chi.FFT("R", "T")
        Factor=self.__Map.Beta/self.__Map.MaxTauBin
        if self.__Lat.Name in ["Square", "Cubic"]:
            StagKIndex=self.__Map.CoordiIndex([e/2 for e in self.__Map.L])
            self.__History["UnifChi"].append(Chi.Data[0,0,0,0,0,0]*Factor)
            self.__History["StagChi"].append(Chi.Data[0,0,0,0,StagKIndex,0]*Factor)
        elif self.__Lat.Name in ["Checkerboard", "3DCheckerboard"]:
            StagKIndex=self.__Map.CoordiIndex([0 for e in self.__Map.L])
            self.__History["UnifChi"].append((Chi.Data[0,0,0,0,0,0]+Chi.Data[0,0,0,1,0,0])*Factor)
            self.__History["StagChi"].append((Chi.Data[0,0,0,0,StagKIndex,0]-Chi.Data[0,0,0,1,StagKIndex,0])*Factor)
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
