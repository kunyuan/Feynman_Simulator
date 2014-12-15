import weight
import lattice as lat
from weight import UP,DOWN,IN,OUT,TAU,SP,SUB,VOL
import numpy as np
from logger import *

class BareFactory:
    def __init__(self, map, Hopping, Interaction, Mu, ExternalField):
        self.__Map=map
        self.__Interaction=Interaction
        self.__ExternalField=ExternalField
        self.__Mu=Mu
        self.__Hopping=Hopping
        self.__MaxTauBin=self.__Map.MaxTauBin
        self.__Beta=self.__Map.Beta
        self.BareG=weight.Weight("G.SmoothT", self.__Map, "TwoSpins", "AntiSymmetric")
        self.BareW=weight.Weight("W.DeltaT", self.__Map, "FourSpins", "Symmetric")

    def Build(self, Model, LatName):
        self.Lat=lat.Lattice(LatName, self.__Map.L)
        if Model=="J1J2" and LatName=="Checkboard":
            self.__J1J2onCheckborad()
        elif Model=="J1J2" and LatName=="Square":
            self.__J1J2onSquare()
        return (self.BareG,self.BareW)

    def __J1J2onCheckborad(self):
        Lx,Ly=self.__Map.L
        #Dimension: 2
        #NSublat: 2
        #Bare G
        self.__Mu=1j*np.pi/2.0/self.__Map.Beta
        self.__Hopping=[0.0]
        log.info("set Mu={0}, Hopping={1}, and SmoothT Bare G".format(self.__Mu, self.__Hopping))
        TauGrid=np.array([self.__Map.IndexToTau(t) for t in range(self.__MaxTauBin)])
        TauWeight=np.exp(self.__Mu*TauGrid)/(1.0+1.0*1j)
        for sp in self.__Map.GetConservedSpinIndexs("TwoSpins"):
            for sub in self.__Map.GetLocalSublatIndexs():
                self.BareG.Data[sp][sub][0][:]=TauWeight[:]
        #Bare W
        J1,J2=self.__Interaction[0:2]
        spinindex=self.__Map.Spin4Index((UP,UP),(UP,UP))
        subA2B=self.__Map.SublatIndex(0,1)
        subB2A=self.__Map.SublatIndex(1,0)
        subA2A=self.__Map.SublatIndex(0,0)
        subB2B=self.__Map.SublatIndex(1,1)
        coordA2B=[(0,0),(0,Ly-1),(Lx-1,0),(Lx-1,Ly-1)]
        coordB2A=[(0, 0),(0,1),(1,0),(1,1)]
        coordA2A=[(0,1),(1,0),(0,Ly-1),(Lx-1,0)]
        coordB2B=[(0, 1),(1,0),(0,Ly-1),(Lx-1,0)]
        #J1 interaction A-->B, B-->A
        for i in coordA2B:
            self.BareW.Data[spinindex,subA2B,self.__Map.CoordiIndex(i)] = J1;
        for i in coordB2A:
            self.BareW.Data[spinindex,subB2A,self.__Map.CoordiIndex(i)] = J1;
        #J2 interaction A-->A, B-->B
        for i in coordA2A:
            self.BareW.Data[spinindex,subA2A,self.__Map.CoordiIndex(i)] = J2;
        for i in coordB2B:
            self.BareW.Data[spinindex,subB2B,self.__Map.CoordiIndex(i)] = J2;
        #Generate other non-zero spin configuration
        for e in self.__Map.GetSpin4SimilarTuples((UP,UP),(UP,UP)):
            self.BareW.Data[self.__Map.Spin4Index(*e),...]=self.BareW.Data[spinindex,...]
        for e in self.__Map.GetSpin4SimilarTuples((DOWN,DOWN),(UP,UP)):
            self.BareW.Data[self.__Map.Spin4Index(*e),...]=-1.0*self.BareW.Data[spinindex,...]
        for e in self.__Map.GetSpin4SimilarTuples((DOWN,UP),(UP,DOWN)):
            self.BareW.Data[self.__Map.Spin4Index(*e),...]=2.0*self.BareW.Data[spinindex,...]

    def __J1J2onSquare(self):
        Lx,Ly=self.__Map.L
        #Dimension: 2
        #NSublat: 2
        #Bare G
        self.__Mu=1j*np.pi/2.0/self.__Map.Beta
        self.__Hopping=[0.0]
        log.info("set Mu={0}, Hopping={1}, and SmoothT Bare G".format(self.__Mu, self.__Hopping))
        TauGrid=np.array([self.__Map.IndexToTau(t) for t in range(self.__MaxTauBin)])
        TauWeight=np.exp(self.__Mu*TauGrid)/(1.0+1.0*1j)
        for sp in self.__Map.GetConservedSpinIndexs("TwoSpins"):
            for sub in self.__Map.GetLocalSublatIndexs():
                self.BareG.Data[sp][sub][0][:]=TauWeight[:]
        #Bare W
        J1,J2=self.__Interaction[0:2]
        spinindex=self.__Map.Spin4Index((UP,UP),(UP,UP))
        sub=self.__Map.SublatIndex(0,0)
        coordnn=[(0,1),(1,0),(Lx-1,0),(0,Ly-1)]
        coordnnn=[(1,1),(Lx-1,1),(1,Ly-1),(Lx-1,Ly-1)]
        #J1 interaction on nearest neighbors
        for i in coordnn:
            self.BareW.Data[spinindex,sub,self.__Map.CoordiIndex(i)] = J1;
        #J2 interaction on next nearest neighbors
        for i in coordnnn:
            self.BareW.Data[spinindex,sub,self.__Map.CoordiIndex(i)] = J2;
        #Generate other non-zero spin configuration
        for e in self.__Map.GetSpin4SimilarTuples((UP,UP),(UP,UP)):
            self.BareW.Data[self.__Map.Spin4Index(*e),...]=self.BareW.Data[spinindex,...]
        for e in self.__Map.GetSpin4SimilarTuples((DOWN,DOWN),(UP,UP)):
            self.BareW.Data[self.__Map.Spin4Index(*e),...]=-1.0*self.BareW.Data[spinindex,...]
        for e in self.__Map.GetSpin4SimilarTuples((DOWN,UP),(UP,DOWN)):
            self.BareW.Data[self.__Map.Spin4Index(*e),...]=2.0*self.BareW.Data[spinindex,...]

    def Plot(self):
        import matplotlib.pyplot as plt
        color=('r','g','b')
        points=self.Lat.GetSitesList(self.__Map)
        for coord, label, sub in points:
            x,y=coord;
            plt.scatter(x,y,s=100,c=color[sub])
            plt.annotate(
                    str(label),
                    xy = (x, y), xytext = (15, 10),
                    textcoords = 'offset points', ha = 'right', va = 'bottom',
                    arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
        lines=self.__GetBareWList()
        for start,end, sub in lines:
            x,y=zip(start,end)
            plt.plot(x,y,c=color[sub],lw=2)
        plt.show()

    def __GetBareWList(self):
        SpinIndex=self.__Map.Spin4Index((UP,UP),(UP,UP))
        offset=np.array(self.__Map.L)/2-1
        BareWList=[]
        for sub in self.__Map.GetAllSublat():
            for coord in self.__Map.GetAllCoordi():
                weight=self.BareW.Data[SpinIndex, self.__Map.SublatIndex(*sub),self.__Map.CoordiIndex(coord)]
                if weight*weight>1.0e-10:
                    BareWList.append([self.Lat.GetRealVec((0,0), sub[IN], offset), \
                                    self.Lat.GetRealVec(coord, sub[OUT], offset), 
                                    sub[IN]])
        return BareWList
