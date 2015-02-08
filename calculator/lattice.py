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
        elif Name=="Kagome":
            self.__Kagome()
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
        self.PathName=["$\Gamma$\n$(0,0)$", "$X$\n(\pi,0)", "$M$\n(\pi,\pi)$","$\Gamma$\n$(0,0)$"]
        self.IndependtBZCenter=[(0,0)]

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
        self.PathName=["$\Gamma$\n$(0,0)$", "$X$\n(\pi,0)", "$M$\n(\pi,\pi)$","$\Gamma$\n$(0,0)$"]
        self.IndependtBZCenter=[(0,0)]

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
        #TODO: self.Path is not correct yet
        self.Path=[(0,0),(PI2/root3,0)]
        self.PathName=["\Gamma", "M"]
    def __Kagome(self):
        self.Dim=2
        self.__AssertDim()
        self.NSublat=3
        self.__AssertNSublat()
        root3=math.sqrt(3)
        PI2=2.0*np.pi
        self.LatVec=np.array([[1.0,0.0],
                             [0.5,root3/2]])
        self.SubLatVec = np.array([[0.5, 0.0],
                                  [1.0/4.0, root3/4],
                                  [3.0/4.0, root3/4]])
        self.ReciprocalLatVec=np.array([[0.0, PI2*2.0/root3],
                                        [PI2, -PI2/root3]])
        #TODO: self.Path is not correct yet
        self.Path=[(0,0),(PI2/root3,0)]
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
        self.Path=[(0,0,0),(PI,0,0),(PI,PI,0),(0,0,0),(PI,PI,PI),(PI,0,0)]
        self.PathName=["$\Gamma$\n$(0,0,0)$", "$X$\n$(\pi,0,0)$", "$M$\n$(\pi,\pi,0)$",
                       "$\Gamma\n$(0,0,0)$", "$R$\n$(\pi,\pi,\pi)$", "$X$\n$(\pi,0,0)$"]
        self.IndependtBZCenter=[(0,0,0)]

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
        self.PathName=["$\Gamma$\n$(0,0,0)$", "$X$\n$(\pi,0,0)$", "$M$\n$(\pi,\pi,0)$",
                       "$\Gamma\n$(0,0,0)$", "$R$\n$(\pi,\pi,\pi)$", "$X$\n$(\pi,0,0)$"]
        self.IndependtBZCenter=[(0,0,0)]

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
        P={"G": (0,0,0), "X":(0,2*PI,0),  "W":(PI,2*PI,0), \
           "K":(1.5*PI,1.5*PI,0),"L": (PI,PI,PI), "U": (PI/2,2*PI,PI/2)}
        L={"G":"$\Gamma$\n$(0,0,0)$", "X":"$X$\n$(0,2\pi,0)$", "W": "$W$\n$(\pi,2\pi,0)$", \
           "K": "$K$\n$(3\pi/2,3\pi/2,0)$", "L": "$L$\n$(\pi,\pi,\pi)$", "U":"$U$\n$(\pi/2,2\pi,0)$"}
        self.Path=[P["G"], P["X"], P["W"], P["K"],
                P["G"], P["L"], P["U"], P["W"], P["L"], P["K"], P["U"], P["X"]]
        self.PathName=[L["G"], L["X"], L["W"], L["K"],
                L["G"], L["L"], L["U"], L["W"], L["L"], L["K"], L["U"], L["X"]]
        self.IndependtBZCenter=[(0,0,0),(2*PI,2*PI,-2*PI),(2*PI,2*PI,2*PI),(4*PI,0,0)]

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

    def FourierTransformation(self, Data, KCoordi, KType="Integer", bound=None):
        """Fourier Transformation in real lattice vector space:
           Data: numpy array with shape [SubLatIn, SubLatOut, VOL]
           KCoordi: list of momentums needed to do fourier transformation, with type KType
           Definition:
               F(K)=1.0/VOL/(Number of Sublattice)*\sum_{i,a;,j,b} G(i,a;j,b)*exp(-i*(r_{i,a}-r_{j,b})*K)
               where i,j are IN/OUT lattice indexies, a,b are IN/OUT sublattice indexies
            """
        DataK=[]
        K=[]
        SiteNum=self.__Map.Vol*self.__Map.NSublat
        vec=np.zeros((SiteNum, self.Dim))
        data=np.zeros(SiteNum)+0*1j
        index=0
        for SubLatIn in range(self.NSublat):
            points,_=self.GetSitesList(HasOffset=False, SubLatIn=SubLatIn)
            for i in range(self.__Map.Vol):
                RealVec, Coord, SubLatOut=points[i]
                vec[index,:]=RealVec
                data[index]=Data[SubLatIn, SubLatOut, self.__Map.CoordiIndex(Coord)]
                index+=1
        Assert(index==SiteNum, "{0} is expected for the total number of sites, not {1}!".format(SiteNum, index))
                
        for p in KCoordi:
            if KType=="Integer":
                KVec=self.ReciprocalLatVec[0,:]*p[0]/self.L[0]
                for i in range(1,self.Dim):
                    KVec+=self.ReciprocalLatVec[i,:]*p[i]/self.L[i]
            elif KType=="Real":
                KVec=np.array(p)
            flag=True
            if bound is not None:
                bound=np.array(bound)
                for d in range(0, bound.shape[0]):
                    if KVec[d]<bound[d][0] or KVec[d]>bound[d][1]:
                        flag=False
                        break
            if flag:
                f=0
                f+=np.dot(data, np.exp(-1j*np.dot(vec[:,:], KVec)))
                K.append(KVec)
                DataK.append(f.real/self.NSublat)
        return K, DataK

    def __Distance(self,a,b):
        return np.sum(np.abs(a-b)**2,axis=-1)**(1./2)
    def __GetKVecBetween(self,start,end,K):
        return K[np.where(np.abs(self.__Distance(K,start)
            +self.__Distance(end,K)-self.__Distance(end,start))<1e-3)]

    def GetKVecAlongPath(self, FromKVec, ToKVec, Center):
        FromKVec=np.array(FromKVec)+np.array(Center)
        ToKVec=np.array(ToKVec)+np.array(Center)
        try:
            flag=np.allclose(self.__Center,Center)
        except:
            flag=False
        finally:
            if not flag:
                #reconstruct KVecList if the BZ center is different from the cached KVecList
                self.__Center=Center
                import itertools 
                KVec=[]
                Index=[]
                for i in range(self.Dim):
                    Index.append(range(-self.L[i], self.L[i]+1))
                KCoordList=list(itertools.product(*Index))
                self.__KVecList=np.zeros((len(KCoordList), self.Dim))+np.array(Center)
                for i in range(len(KCoordList)):
                    for j in range(self.Dim):
                        self.__KVecList[i,:] += KCoordList[i][j]*self.ReciprocalLatVec[j,:]/self.L[j]
        KList=list(self.__GetKVecBetween(FromKVec, ToKVec, self.__KVecList))
        return KList

    def GetRealVec(self, Coordi, SubLatIn, SubLatOut, offset):
        '''
           Coordi: D-dimensional vector of coordinates
           SubLat: only the OUT sublattice is needed, IN sublattice is assumed to be 0
        '''
        v=self.__Shift(Coordi+offset)
        return tuple(np.einsum("ij,i->j",self.LatVec,v)+self.SubLatVec[SubLatOut]-self.SubLatVec[SubLatIn])
    def GetSitesList(self, HasOffset=True, SubLatIn=0):
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
            for SubLatOut in range(self.NSublat):
                Points[self.__Map.LatIndex(coord, SubLatOut)]=[self.GetRealVec(coord,SubLatIn,SubLatOut,offset), \
                        coord,SubLatOut]
                for subN in range(SubLatOut, self.NSublat):
                    LinesInUnitCell.append([(self.__Map.LatIndex(coord, SubLatOut), \
                                              self.__Map.LatIndex(coord, subN)), SubLatOut])
        return Points, LinesInUnitCell

if __name__=="__main__":
    import weight
    WeightPara={"NSublat": 1, "L":[4, 4, 4],
            "Beta": 0.5, "MaxTauBin":64}
    Map=weight.IndexMap(**WeightPara)
    l=Lattice("Cubic", Map)

    with open("Coordinates.txt","w") as f:
        f.write(str(l.GetSitesList()))
