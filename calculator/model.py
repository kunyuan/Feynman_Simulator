import weight
import lattice as lat
from weight import UP,DOWN,IN,OUT
import numpy as np
from logger import *

class BareFactory:
    def __init__(self, map, Hamiltonian):
        self.__Map=map
        self.__Interaction=np.array(Hamiltonian["Interaction"])
        self.__ExternalField=np.array(Hamiltonian["ExternalField"])
        self.__Mu=np.array(Hamiltonian["ChemicalPotential"])
        self.__Hopping=np.array(Hamiltonian["Hopping"])
        self.__MaxTauBin=self.__Map.MaxTauBin
        self.__Beta=self.__Map.Beta
        self.BareG=weight.Weight("SmoothT", self.__Map, "TwoSpins", "AntiSymmetric")
        self.BareW=weight.Weight("DeltaT", self.__Map, "FourSpins", "Symmetric")

    def Build(self, Model, LatName=None):
        self.Lat=lat.Lattice(LatName, self.__Map)
        if Model=="Heisenberg":
            self.__Heisenberg(LatName)
        if Model=="J1J2":
            self.__J1J2(LatName)
        elif Model=="DiagCount":
            self.__DiagCount()
        return (self.BareG,self.BareW)

    def __Heisenberg(self, LatName):
        Assert(len(self.__Interaction)==1, "Heisenberg model only has one coupling!")
        self.__SpinModel(LatName)
    def __J1J2(self, LatName):
        Assert(len(self.__Interaction)==2, "J1J2 model only has two couplings!")
        self.__SpinModel(LatName)

    def __SpinModel(self, LatName):
        Beta=self.__Map.Beta
        #Bare G
        self.__Hopping=np.array([0.0])
        log.info("set Mu={0}, Hopping={1}, and SmoothT Bare G".format(self.__Mu, self.__Hopping))
        Assert(len(self.__ExternalField)>=self.__Map.NSublat, 
                "expect at least {0} externalfield components!".format(self.__Map.NSublat))

        TauGrid=np.array([self.__Map.IndexToTau(t) for t in range(self.__MaxTauBin)])
        Pauli_Z=self.__Map.Pauli()[2]
        for sub in range(self.__Map.NSublat):
            for sp in range(2):
                Mu=1j*np.pi/2.0/Beta+Pauli_Z[sp, sp]*self.__ExternalField[sub]
                self.BareG.Data[sp,sub,sp,sub,0,:]=np.exp(Mu*TauGrid)/(1.0+np.exp(Mu*Beta))

        Interaction=self.__Interaction+[0,0,0,0,0]
        J1,J2=Interaction[0:2]
        #Bare W
        #Dimension: 2
        spin=self.__Map.Spin2Index(UP,UP)
        if LatName=="Checkboard":
        #NSublat: 2
            Lx,Ly=self.__Map.L
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
            Lx,Ly=self.__Map.L
            sub=0
            coordnn=[(0,1),(1,0),(Lx-1,0),(0,Ly-1)]
            coordnnn=[(1,1),(Lx-1,1),(1,Ly-1),(Lx-1,Ly-1)]
            #J1 interaction on nearest neighbors
            for i in coordnn:
                self.BareW.Data[spin,sub,spin,sub,self.__Map.CoordiIndex(i)] = J1/4.0;
            #J2 interaction on next nearest neighbors
            for i in coordnnn:
                self.BareW.Data[spin,sub,spin,sub,self.__Map.CoordiIndex(i)] = J2/4.0;
        elif LatName=="Cubic":
        #NSublat: 1
            Lx,Ly,Lz=self.__Map.L
            sub=0
            coordnn=[(1,0,0),(0,1,0),(0,0,1),(Lx-1,0,0),(0,Ly-1,0),(0,0,Lz-1)]
            #coordnnn=[(1,1),(Lx-1,1),(1,Ly-1),(Lx-1,Ly-1)]
            coordnnn=[]
            #J1 interaction on nearest neighbors
            for i in coordnn:
                self.BareW.Data[spin,sub,spin,sub,self.__Map.CoordiIndex(i)] = J1/4.0;
            #J2 interaction on next nearest neighbors
            for i in coordnnn:
                self.BareW.Data[spin,sub,spin,sub,self.__Map.CoordiIndex(i)] = J2/4.0;
            pass
        if LatName=="Pyrochlore":
        #NSublat: 4
            Lx,Ly,Lz=self.__Map.L
            A,B,C,D=0,1,2,3
            coord=[]
            for i in range(4):
                coord.append([])
                for j in range(4):
                    coord[i].append([])
            coord[A][B]=[(0,0,0),(0,Ly-1,0)]
            coord[A][C]=[(0,0,0),(0,0,Lz-1)]
            coord[A][D]=[(0,0,0),(Lx-1,0,0)]
            coord[B][A]=[(0,0,0),(0,1,0)]
            coord[B][C]=[(0,0,0),(0,1,Lz-1)]
            coord[B][D]=[(0,0,0),(Lx-1,1,0)]
            coord[C][A]=[(0,0,0),(0,0,1)]
            coord[C][B]=[(0,0,0),(0,Ly-1,1)]
            coord[C][D]=[(0,0,0),(Lx-1,0,1)]
            coord[D][A]=[(0,0,0),(1,0,0)]
            coord[D][B]=[(0,0,0),(1,Ly-1,0)]
            coord[D][C]=[(0,0,0),(1,0,Lz-1)]

            for i in range(4):
                for j in range(4):
                    for e in coord[i][j]:
                        self.BareW.Data[spin,i,spin,j,self.__Map.CoordiIndex(e)] = J1/4;
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

    def ToDict(self):
        points, LinesInUnitCell=self.Lat.GetSitesList()
        interaction=self.__GetBareWList()
        return {"Points":points, "Interaction":interaction, "Lines": LinesInUnitCell}

    def Plot(self):
        if self.Lat.Dim==2:
            import matplotlib.pyplot as plt
            Assert(self.Lat.Dim==2, "Plot() only works for two dimensional system for now.")
            color=('r','g','b')
            points, _=self.Lat.GetSitesList()
            for vec, coord, sub in points:
                x,y=vec;
                plt.scatter(x,y,s=100,c=color[sub])
                plt.annotate(
                        str(coord),
                        xy = (x, y), xytext = (15, 10),
                        textcoords = 'offset points', ha = 'right', va = 'bottom',
                        arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
            for vec, coord, sub in self.__GetBareWList():
                start,end=vec
                x,y=zip(start, end)
                plt.plot(x,y,c=color[sub],lw=2)
            plt.show()

    def __GetBareWList(self):
        if self.Lat.Dim==2:
            Origin=(0,0)
        elif self.Lat.Dim==3:
            Origin=(0,0,0)
        Spin=self.__Map.Spin2Index(UP,UP)
        offset=np.array(self.__Map.L)/2-1
        size=self.__Map.Vol*self.__Map.NSublat
        BareWList=[]
        for sub in self.__Map.GetAllSublatTuple():
            for coord in self.__Map.GetAllCoordi():
                weight=self.BareW.Data[Spin,sub[IN],Spin,sub[OUT],self.__Map.CoordiIndex(coord)]
                if weight*weight>1.0e-10:
                    n=self.__Map.LatIndex(coord, sub[OUT])
                    vec=(self.Lat.GetRealVec(Origin, sub[IN], offset), \
                                    self.Lat.GetRealVec(coord, sub[OUT], offset))
                    coord=(self.__Map.LatIndex(Origin, sub[IN]), \
                                      self.__Map.LatIndex(coord, sub[OUT]))
                    BareWList.append([vec, coord, sub[IN]])
        return BareWList
