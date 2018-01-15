#!/usr/bin/env python
import numpy as np
import numpy.linalg as linalg
import IO
# from logger import *
# import parameter as para

class SVDBasis:
    def __init__(self, MaxTauBin, beta, bosefermi):
        self.__MaxTauBin=MaxTauBin
        self.__Beta=beta
        self.__Basis={"Beta":beta,"MaxTauBin":MaxTauBin, "BoseFermi":bosefermi}

    def FermiKernel(self, w, t, beta):
        x=beta*w/2
        y=2*t/beta-1
        if x>100:
            return np.exp(-x*(y+1.))
        if x<-100:
            return np.exp(x*(1.0-y))
        return np.exp(-x*y)/(2*np.cosh(x))

    def BoseKernel(self, w, t, beta):
        x=beta*w/2
        y=2*t/beta-1
        if x>200:
            return w*np.exp(-x*(y+1.))
        if x<-200:
            return -w*np.exp(x*(1.0-y))
        return w*np.exp(-x*y)/(2*np.sinh(x))

    def GenerateBasis(self, N):
        Nw=1000
        w=np.linspace(-100,100,Nw)
        Nt=self.__MaxTauBin
        t=np.linspace(0, self.__Beta, Nt+1)
        t=np.array([e+1.0/Nt/2 for e in t[:-1]])
        kMatrix=np.zeros([Nw,Nt])
        if self.__Basis["BoseFermi"] == "Fermi":
            for i in range(len(w)):
                kMatrix[i,:]=self.FermiKernel(w[i],t,self.__Beta)
        elif self.__Basis["BoseFermi"] == "Bose":
            for i in range(len(w)):
                kMatrix[i,:]=self.BoseKernel(w[i],t,self.__Beta)

        u,s,v=linalg.svd(kMatrix)
        v_inv=linalg.inv(v)

        self.__Basis["Number"]=N
        self.__Basis["Basis"]=[]
        for i in range(N):
            self.__Basis["Basis"].append(list(v_inv[:,i]))
    
    def getBasis(self):
        return self.__Basis

    def Save(self, filename, mode="w"):
        IO.SaveDict(filename, mode, self.__Basis)

if __name__=="__main__":
    beta=0.8
    MaxTauBin=128
    svd=SVDBasis(MaxTauBin,beta, "Fermi")
    svd.GenerateBasis(30)
    svd.Save("FermiBasis.dat")
    svd=SVDBasis(MaxTauBin,beta, "Bose")
    svd.GenerateBasis(30)
    svd.Save("BoseBasis.dat")












