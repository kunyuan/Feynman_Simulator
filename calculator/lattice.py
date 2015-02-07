#!usr/bin/env python
import numpy as np
import math
from logger import *
PI=np.pi

class Lattice:
    def __init__(self, Name, Map):
        self.__Map=Map
        self.L=np.array(Map.L)
        self.Name=Name
        if Name=="Checkerboard":
            self.__Checkerboard()
        elif Name=="3DCheckerboard":
            self.__3DCheckerboard()
        elif Name=="Honeycomb":
            self.__Honeycomb()
        elif Name=="Square":
            self.__Square()
        elif Name=="Cubic":
            self.__Cubic()
        elif Name=="Pyrochlore":
            self.__Pyrochlore()
        else:
            Assert(False, "Lattice {0} has not been implemented!".format(self.Name))

    #2D lattice
    def __Checkerboard(self):
        self.Dim=2
        self.__AssertDim()
        self.NSublat=2
        self.__AssertNSublat()
        self.LatVec=np.array([[1.0,1.0],
                              [1.0,-1.0]])
        self.SubLatVec=np.array([[0.0,0.0],
                                 [1.0,0.0]])
        self.ReciprocalLatVec =np.array([[PI, -PI],
                                         [PI,  PI]])
        self.Path=[(0,0),(PI,0),(PI,PI),(0,0)]
        self.PathNum=[self.L[0]/2, self.L[1]/2, self.L[0]/2]
        self.PathName=["\Gamma", "\Chi", "M","\Gamma"]

    def __Square(self):
        self.Dim=2
        self.__AssertDim()
        self.NSublat=1
        self.__AssertNSublat()
        self.LatVec=np.array([[1.0,0.0],
                              [0.0,1.0]])
        self.SubLatVec=np.array([[0.0,0.0]])
        self.ReciprocalLatVec =np.array([[2.0 * PI, 0.0],
                                         [0.0, 2.0 * PI]])
        self.Path=[(0,0),(PI,0),(PI,PI),(0,0)]
        self.PathNum=[self.L[0]/2, self.L[1]/2, self.L[0]/2]
        self.PathName=["\Gamma", "\Chi", "M","\Gamma"]

    def __Honeycomb(self):
        self.Dim=2
        self.__AssertDim()
        self.NSublat=2
        self.__AssertNSublat()
        root3=math.sqrt(3)
        PI2=2.0*np.pi
        self.LatVec=np.array([[0.0,1.0],
                             [root3/2,-0.5]])
        self.SubLatVec = np.array([[0.0, 0.0],
                                  [1.0/2.0/root3, 0.5]])
        self.ReciprocalLatVec=np.array([[PI2/root3,PI2],
                                       [PI2*2.0/root3,0.0]])
        self.Path=[(0,0),(PI2/root3,0)]
        self.PathNum=[self.L[0]/2]
        self.PathName=["\Gamma", "M"]
    #3D lattice
    def __Cubic(self):
        self.Dim=3
        self.__AssertDim()
        self.NSublat=1
        self.__AssertNSublat()
        self.LatVec=np.array([[1.0,0.0,0.0],
                              [0.0,1.0,0.0],
                              [0.0,0.0,1.0]])
        self.SubLatVec=np.array([[0.0,0.0,0.0]])
        self.ReciprocalLatVec =np.array([[2.0 * np.pi, 0.0, 0.0],
                                         [0.0, 2.0 * np.pi, 0.0],
                                         [0.0, 0.0, 2.0 * np.pi]])

    def __3DCheckerboard(self):
        self.Dim=3
        self.__AssertDim()
        self.NSublat=2
        self.__AssertNSublat()
        self.LatVec=np.array([[1.0,1.0,0.0],
                              [1.0,-1.0,0.0],
                              [1.0,0.0,1.0]])
        self.SubLatVec=np.array([[0.0,0.0,0.0],
                                 [1.0,0.0,0.0]])
        self.ReciprocalLatVec =np.array([[PI, PI,-PI],
                                         [PI,-PI,-PI],
                                         [ 0, 0,2*PI]])
        self.Path=[(0,0,0),(PI,0,0),(PI,PI,0),(0,0,0),(PI,PI,PI),(PI,0,0)]
        self.PathNum=[self.L[0]/2, self.L[1]/2, self.L[0]/2, self.L[0]/2, self.L[2]/2]
        self.PathName=["\Gamma", "\Chi", "M","\Gamma", "\R", "\Chi"]

    def __Pyrochlore(self):
        self.Dim=3
        self.__AssertDim()
        self.NSublat=4
        self.__AssertNSublat()
        self.LatVec=np.array([[0.5,0.0,0.5],
                              [0.5,0.5,0.0],
                              [0.0,0.5,0.5]])
        self.SubLatVec=np.array([[0.0,0.0,0.0],
                                 [0.0,0.25,0.25],
                                 [0.25,0.0,0.25],
                                 [0.25,0.25,0]])
        self.ReciprocalLatVec=np.array([[2*PI, -2*PI, 2*PI],
                                        [2*PI, 2*PI, -2*PI],
                                        [-2*PI, 2*PI, 2*PI]])
        self.Path=[(0,2*PI,0),(PI/2,2*PI,PI/2),(PI,PI,PI),(0,0,0),(0,2*PI,0),(PI,2*PI,0), \
                (3*PI/2,3*PI/2,0)]
        self.PathNum=[self.L[0]/8, self.L[0]/8, self.L[0]/2, self.L[0]/2, self.L[0]/4, self.L[1]/8]
        self.PathName=["X", "U", "L","\Gamma", "X", "W", "K"]


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

    def FourierTransformation_RealSpace(self, Data, KCoordi, KType="Integer"):
        DataK=[]
        K=[]
        points,_=self.GetSitesList(False)
        #for vec, coord, sub in points:
            #if sub==1:
                #LatPoints.append((np.array(vec), self.__Map.CoordiIndex(coord), sub))
        vec=np.zeros((len(points), self.Dim))
        data=np.zeros(len(points))+0*1j
        for i in range(len(points)):
        #for vec, coord, sub in points:
            if points[i][2]==0:
                vec[i,:]=points[i][0]
                data[i]=Data[points[i][2], self.__Map.CoordiIndex(points[i][1])]
                
        for p in KCoordi:
            if KType=="Integer":
                KVec=self.ReciprocalLatVec[0,:]*p[0]/self.L[0]
                for i in range(1,self.Dim):
                    KVec+=self.ReciprocalLatVec[i,:]*p[i]/self.L[i]
            elif KType=="Real":
                KVec=np.array(p)
            f=0
            f+=np.dot(data, np.exp(-1j*np.dot(vec[:,:], KVec)))
            #for vec, coord, sub in LatPoints:
                #if sub==1:
                    #f+=Data[sub,coord]*np.exp(-1j*np.dot(vec,KVec))
            K.append(KVec)
            DataK.append(f.real)
        return K, DataK

    def GetRealVec(self, Coordi, SubLat, offset):
        '''
           Coordi: D-dimensional vector of coordinates
           SubLat: only the OUT sublattice is needed, IN sublattice is assumed to be 0
        '''
        v=self.__Shift(Coordi+offset)
        return tuple(np.einsum("ij,i->j",self.LatVec,v)+self.SubLatVec[SubLat])
    def GetSitesList(self, HasOffset=True):
        """
        return: list of all sites, with format 
                [tuple of real space vectors of sites, tuple of integer vectors of coordinates, SubLat] 
        """
        if HasOffset:
            offset=self.L/2-1
        else:
            offset=np.array([0 for e in self.L])
        Points=[None,]*(self.__Map.Vol*self.__Map.NSublat)
        LinesInUnitCell=[]
        #Origin=[0 for e in self.L]
        for coord in self.__Map.GetAllCoordi():
            for sub in range(self.NSublat):
                Points[self.__Map.LatIndex(coord, sub)]=[tuple(self.GetRealVec(coord,sub,offset)),coord,sub]
                for subN in range(sub, self.NSublat):
                    LinesInUnitCell.append([(self.__Map.LatIndex(coord, sub), \
                                              self.__Map.LatIndex(coord, subN)), sub])
        return Points, LinesInUnitCell

if __name__=="__main__":
    import weight
    WeightPara={"NSublat": 1, "L":[4, 4, 4],
            "Beta": 0.5, "MaxTauBin":64}
    Map=weight.IndexMap(**WeightPara)
    l=Lattice("Cubic", Map)

    with open("Coordinates.txt","w") as f:
        f.write(str(l.GetSitesList()))
