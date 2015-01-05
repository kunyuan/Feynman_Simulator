#!usr/bin/env python
import numpy as np
import math
from logger import *

class Lattice:
    def __init__(self, Name, Map):
        self.__Map=Map
        self.L=np.array(Map.L)
        self.Name=Name
        print Name
        if Name=="Checkboard":
            self.__Checkboard()
        elif Name=="Honeycomb":
            self.__Honeycomb()
        elif Name=="Square":
            self.__Square()
        else:
            Assert(False, "Not implemented!")

    def __Checkboard(self):
        self.Dim=2
        self.__AssertDim()
        self.NSublat=2
        self.__AssertNSublat()
        self.LatVec=np.array([[1.0,0.0],
                              [0.0,1.0]])
        self.SubLatVec=np.array([[0.0,0.0],
                                 [0.5,0.5]])
        self.ReciprocalLatVec =np.array([[2.0 * np.pi, 0.0],
                                         [0.0, 2.0 * np.pi]])
    def __Square(self):
        self.Dim=2
        self.__AssertDim()
        self.NSublat=1
        self.__AssertNSublat()
        self.LatVec=np.array([[1.0,0.0],
                              [0.0,1.0]])
        self.SubLatVec=np.array([[0.0,0.0]])
        self.ReciprocalLatVec =np.array([[2.0 * np.pi, 0.0],
                                         [0.0, 2.0 * np.pi]])
    def __Honeycomb(self):
        self.Dim=2
        self.__AssertDim()
        self.NSublat=2
        self.__AssertNSublat()
        root3=math.sqrt(3)
        PI2=2.0*np.pi
        self.LatVec=([[0.0,1.0],
                      [root3/2,-0.5]])
        self.SubLatVec=([[PI2/root3,PI2],
                         [PI2*2.0/root3,0.0]])
        self.ReciprocalLatVec = ([[0.0, 0.0],
                                  [1.0/2.0/root3, 0.5]])

    def __AssertDim(self):
        Assert(len(self.L)==self.Dim, "Dimension {0} is expected for {1} Lattice, not {2}" \
                .format(self.Dim, self.Name, len(self.L)))
    def __AssertNSublat(self):
        Assert(self.NSublat==self.__Map.NSublat, "{0} is expected for {1} lattice, not {2}"
                .format(self.NSublat, self.Name, self.__Map.NSublat))
    def __Shift(self,Coordi):
        v=list(Coordi) #make a copy of Vec
        for i in range(len(Coordi)):
            if v[i]<0.0:
                v[i]+=self.L[i]
            if v[i]>=self.L[i]:
                v[i]-=self.L[i]
        return v
    def GetRealVec(self, Coordi, SubLat, offset):
        '''
           Coordi: D-dimensional vector of coordinates
           SubLat: only the OUT sublattice is needed, IN sublattice is assumed to be 0
        '''
        v=self.__Shift(Coordi+offset)
        return np.einsum("ij,i->j",self.LatVec,v)+self.SubLatVec[SubLat]
    def GetSitesList(self):
        """
        return: list of all sites, with format 
                [tuple of real space vectors of sites, tuple of integer vectors of coordinates, SubLat] 
        """
        offset=self.L/2-1
        Points=[]
        for sub in range(self.NSublat):
            for coord in self.__Map.GetAllCoordi():
                Points.append([tuple(self.GetRealVec(coord,sub,offset)),coord,sub])
        return Points

if __name__=="__main__":
    l=Lattice("Checkboard", [4,4])
    print l.GetSitesList()
