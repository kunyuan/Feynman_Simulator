import weight
import lattice as lat
from weight import UP,DOWN,IN,OUT,TAU,SP1,SUB1,SP2,SUB2,VOL
import numpy as np
from logger import *

class BareFactory:
    def __init__(self, map, Hamiltonian):
        self.__Map=map
        self.__Interaction=Hamiltonian["Interaction"]
        self.__ExternalField=Hamiltonian["ExternalField"]
        self.__Mu=Hamiltonian["ChemicalPotential"]
        self.__Hopping=Hamiltonian["Hopping"]
        self.__MaxTauBin=self.__Map.MaxTauBin
        self.__Beta=self.__Map.Beta
        self.BareG=weight.Weight("SmoothT", self.__Map, "TwoSpins", "AntiSymmetric")
        self.BareW=weight.Weight("DeltaT", self.__Map, "FourSpins", "Symmetric")

    def Build(self, Model, LatName=None):
        self.Lat=lat.Lattice(LatName, self.__Map)
        if Model=="J1J2":
            self.__J1J2(LatName)
        elif Model=="DiagCount":
            self.__DiagCount()
        return (self.BareG,self.BareW)

    def __J1J2(self, LatName):
        Lx,Ly=self.__Map.L
        #Dimension: 2
        #Bare G
        self.__Mu=1j*np.pi/2.0/self.__Map.Beta
        self.__Hopping=[0.0]
        log.info("set Mu={0}, Hopping={1}, and SmoothT Bare G".format(self.__Mu, self.__Hopping))
        TauGrid=np.array([self.__Map.IndexToTau(t) for t in range(self.__MaxTauBin)])
        TauWeight=np.exp(self.__Mu*TauGrid)/(1.0+1.0*1j)
        for sp in self.__Map.GetConservedSpinTuple("TwoSpins"):
            for sub in self.__Map.GetLocalSublatTuple():
                self.BareG.Data[sp[IN]][sub[IN]][sp[OUT]][sub[OUT]][0][:]=TauWeight[:]

        #Bare W
        J1,J2=self.__Interaction[0:2]
        spin=self.__Map.Spin2Index(UP,UP)
        if LatName=="Checkboard":
        #NSublat: 2
            subA=0
            subB=1
            coordA2B=[(0, 0),(0,Ly-1),(Lx-1,0),(Lx-1,Ly-1)]
            coordB2A=[(0, 0),(0,   1),(1,   0),(   1,   1)]
            coordA2A=[(0, 1),(1,   0),(0,Ly-1),(Lx-1,   0)]
            coordB2B=[(0, 1),(1,   0),(0,Ly-1),(Lx-1,   0)]

            #J1 interaction A-->B, B-->A
            for i in coordA2B:
                self.BareW.Data[spin,subA,spin,subB,self.__Map.CoordiIndex(i)] = J1/4;

            for i in coordB2A:
                self.BareW.Data[spin,subB,spin,subA,self.__Map.CoordiIndex(i)] = J1/4;

            #J2 interaction A-->A, B-->B
            for i in coordA2A:
                self.BareW.Data[spin,subA,spin,subA,self.__Map.CoordiIndex(i)] = J2/4;
            for i in coordB2B:
                self.BareW.Data[spin,subB,spin,subB,self.__Map.CoordiIndex(i)] = J2/4;
        elif LatName=="Square":
        #NSublat: 1
            sub=0
            coordnn=[(0,1),(1,0),(Lx-1,0),(0,Ly-1)]
            coordnnn=[(1,1),(Lx-1,1),(1,Ly-1),(Lx-1,Ly-1)]
            #J1 interaction on nearest neighbors
            for i in coordnn:
                self.BareW.Data[spin,sub,spin,sub,self.__Map.CoordiIndex(i)] = J1/4.0;
            #J2 interaction on next nearest neighbors
            for i in coordnnn:
                self.BareW.Data[spin,sub,spin,sub,self.__Map.CoordiIndex(i)] = J2/4.0;
        else:
            Assert(False, "Not implemented yet!")

        #Generate other non-zero spin configuration
        for e in self.__Map.GetSpin4SimilarTuples((UP,UP),(UP,UP)):
            spleft = self.__Map.Spin2Index(*e[IN]) 
            spright = self.__Map.Spin2Index(*e[OUT]) 
            self.BareW.Data[spleft,:,spright,...]=self.BareW.Data[spin,:,spin,...]
        for e in self.__Map.GetSpin4SimilarTuples((DOWN,DOWN),(UP,UP)):
            spleft = self.__Map.Spin2Index(*e[IN]) 
            spright = self.__Map.Spin2Index(*e[OUT]) 
            self.BareW.Data[spleft,:,spright,...]=-1.0*self.BareW.Data[spin,:,spin,...]
        for e in self.__Map.GetSpin4SimilarTuples((DOWN,UP),(UP,DOWN)):
            spleft = self.__Map.Spin2Index(*e[IN]) 
            spright = self.__Map.Spin2Index(*e[OUT]) 
            self.BareW.Data[spleft,:,spright,...]=2.0*self.BareW.Data[spin,:,spin,...]

    def Plot(self):
        import matplotlib.pyplot as plt
        Assert(self.Lat.Dim==2, "Plot() only works for two dimensional system for now.")
        color=('r','g','b')
        points=self.Lat.GetSitesList()
        for vec, coord, sub in points:
            x,y=vec;
            plt.scatter(x,y,s=100,c=color[sub])
            plt.annotate(
                    str(coord),
                    xy = (x, y), xytext = (15, 10),
                    textcoords = 'offset points', ha = 'right', va = 'bottom',
                    arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
        lines=self.__GetBareWList()
        for start,end, sub in lines:
            x,y=zip(start,end)
            plt.plot(x,y,c=color[sub],lw=2)
        plt.show()

    def __GetBareWList(self):
        Spin=self.__Map.Spin2Index(UP,UP)
        offset=np.array(self.__Map.L)/2-1
        BareWList=[]
        for sub in self.__Map.GetAllSublatTuple():
            for coord in self.__Map.GetAllCoordi():
                weight=self.BareW.Data[Spin,sub[IN],Spin,sub[OUT],self.__Map.CoordiIndex(coord)]
                if weight*weight>1.0e-10:
                    BareWList.append([self.Lat.GetRealVec((0,0), sub[IN], offset), \
                                    self.Lat.GetRealVec(coord, sub[OUT], offset), 
                                    sub[IN]])
        return BareWList
