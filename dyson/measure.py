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

    def Measure(self, Chi, BKChi, Determinate, G, NN):
        self.Append("1-JP", Determinate.min())
        Factor=self.__Map.Beta/self.__Map.MaxTauBin

        Chi.FFT("R", "W")
        BKChi.FFT("R", "W")
        UnifK=(0.0,)*self.__Map.Dim
        _, ChiK=self.__Lat.FourierTransformation(Chi.Data[0,:,0,:,:,0]*Factor, [UnifK,],"Real")
        _, BKChiK=self.__Lat.FourierTransformation(BKChi.Data[0,:,0,:,:,0]*Factor, [UnifK,],"Real")
        # print "BKChiK"
        # print BKChiK.shape
        self.Append("UnifChi", ChiK[0])
        self.Append("BKUnifChi", BKChiK[0])
        self.Append("Mag^2 density", ChiK[0]/self.__Map.Vol/self.__Map.NSublat/self.__Map.Beta)
        self.Append("BK Mag^2 density", BKChiK[0]/self.__Map.Vol/self.__Map.NSublat/self.__Map.Beta)

        Chi.FFT("R","T")
        BKChi.FFT("R","T")
        energy=0j
        BK_energy=0j
        for i in range(self.__Lat.NSublat):
            for j in range(self.__Lat.NSublat):
                for l in NN[i][j]:
                    energy+=Chi.Data[0,i,0,j,self.__Map.CoordiIndex(l),0]/self.__Lat.NSublat
                    BK_energy+=BKChi.Data[0,i,0,j,self.__Map.CoordiIndex(l),0]/self.__Lat.NSublat
        self.Append("Energy", energy)
        self.Append("SumRule", Chi.Data[0,0,0,0,0,0])

        self.Append("BK_Energy", BK_energy)
        self.Append("BK_SumRule", BKChi.Data[0,0,0,0,0,0])


        if self.__Lat.Name in ["Chain", "Assymetric_Triangular", "Square", "Cubic", "Checkerboard", "3DCheckerboard", "ValenceBond"]:
            # Chi.FFT("R", "T")
            # print Chi.Data[0,0,0,0,:,0]
            StagK=(PI,)*self.__Map.Dim
            _, ChiK=self.__Lat.FourierTransformation(Chi.Data[0,:,0,:,:,0]*Factor, [StagK,],"Real")
            _, BKChiK=self.__Lat.FourierTransformation(BKChi.Data[0,:,0,:,:,0]*Factor, [StagK,],"Real")
            self.Append("StagChi", ChiK[0])
            self.Append("BKStagChi", BKChiK[0])
            self.Append("StagMag^2 density", ChiK[0]/self.__Map.Vol/self.__Map.NSublat/self.__Map.Beta)
            self.Append("BK StagMag^2 density", BKChiK[0]/self.__Map.Vol/self.__Map.NSublat/self.__Map.Beta)

            #Chi.FFT("R", "T")
            #_, ChiK=self.__Lat.FourierTransformation(Chi.Data[0,:,0,:,:,0], [UnifK,StagK],"Real")
            #self.Append("Mag^2 density", ChiK[0]/self.__Map.Vol/self.__Map.NSublat)
            #self.Append("StagMag^2 density", ChiK[1]/self.__Map.Vol/self.__Map.NSublat)

        elif self.__Lat.Name in ["Triangular"]:
            Chi.FFT("R", "W")
            BKChi.FFT("R", "W")
            offset=0
            x=[]
            KList=[]
            BZcenter=(0,0)
            start, end=np.array(self.__Lat.Path[0]),np.array(self.__Lat.Path[1])
            for k in self.__Lat.GetKVecAlongPath(start, end, BZcenter):
                pos=offset+np.linalg.norm(k-np.array(BZcenter)-start)
                x.append(pos)
                KList.append(k)
            offset+=np.linalg.norm(end-start)
            StagK=KList[-1]
            _, ChiK=self.__Lat.FourierTransformation(Chi.Data[0,:,0,:,:,0]*Factor, [StagK,],"Real")
            _, BKChiK=self.__Lat.FourierTransformation(BKChi.Data[0,:,0,:,:,0]*Factor, [StagK,],"Real")
            self.Append("StagChi", ChiK[0])
            self.Append("BKStagChi", BKChiK[0])
            self.Append("StagMag^2 density", ChiK[0]/self.__Map.Vol/self.__Map.NSublat/self.__Map.Beta)
            self.Append("BK StagMag^2 density", BKChiK[0]/self.__Map.Vol/self.__Map.NSublat/self.__Map.Beta)

        elif self.__Lat.Name in ["Pyrochlore"]:
            Chi.FFT("R", "W")
            KList=[(0.0,0.0,0.0), (4*PI, 2*PI ,0)] #High symmetry point with strongest order
            _, ChiK=self.__Lat.FourierTransformation(Chi.Data[0,:,0,:,:,0]*Factor, KList,"Real")
            self.Append("UnifChi", ChiK[0])
            self.Append("Chi_X(4Pi,2Pi,0)", ChiK[1])

        elif self.__Lat.Name in ["DecoPyrochlore"]:
            Chi.FFT("R", "W")
            KList=[(0.0,0.0,0.0), (4*PI, 2*PI ,0)] #High symmetry point with strongest order
            _, ChiK=self.__Lat.FourierTransformation(Chi.Data[0,:,0,:,:,0]*Factor, KList,"Real")
            self.Append("UnifChi", ChiK[0])
            self.Append("Chi_X(4Pi,2Pi,0)", ChiK[1])


        elif self.__Lat.Name in ["Kagome", "Honeycomb","Triangular"]:
            Chi.FFT("R", "W")
            # KList=[(0.0,0.0),] #High symmetry point with strongest order
            # _, ChiK=self.__Lat.FourierTransformation(Chi.Data[0,:,0,:,:,0]*Factor, KList,"Real")
            # self.Append("UnifChi", ChiK[0])

        else:
            Assert(False, "model not implemented!")

        if self.__Lat.Name in ["ValenceBond"]:
            Chi.FFT("R","W")
            Neighbor1=self.__Map.CoordiIndex((0,0))
            Neighbor2=self.__Map.CoordiIndex((0,self.__Lat.L[1]-1))
            self.Append("VBS", (Chi.Data[0,0,0,1,Neighbor1,0]-Chi.Data[0,0,0,1,Neighbor2,0])*Factor)

        G.FFT("R","T")
        for i in range(self.__Map.NSublat):
            self.Append("<Sz_{0}>".format(i), -0.5*(1.5*G.Data[UP,i,UP,i,0,-1]-0.5*G.Data[UP,i,UP,i,0,-2]-1.5*G.Data[DOWN,i,DOWN,i,0,-1]+0.5*G.Data[DOWN,i,DOWN,i,0,-2]))

        G.FFT("K","T")
        for i in range(self.__Map.NSublat):
            self.Append("<n_{0}>".format(i), np.sum(1.5*G.Data[UP,i,UP,i,:,-1]-0.5*G.Data[UP,i,UP,i,:,-2]+1.5*G.Data[DOWN,i,DOWN,i,:,-1]-0.5*G.Data[DOWN,i,DOWN,i,:,-2])/self.__Map.Vol)

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
