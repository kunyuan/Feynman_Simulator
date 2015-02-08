import weight
import lattice as lat
from weight import UP,DOWN,IN,OUT
import numpy as np
import math
from logger import *

class BareFactory:
    def __init__(self, map, Lat, Hamiltonian, Anneal):
        self.Lat=Lat
        self.__Map=map
        self.__Model=Hamiltonian["Name"]
        self.__Interaction=np.array(Hamiltonian["Interaction"])
        self.__ExternalField=np.array(Hamiltonian["ExternalField"])
        self.__DeltaField=np.array(Anneal["DeltaField"])
        self.__Mu=np.array(Hamiltonian["ChemicalPotential"])
        self.__Hopping=np.array(Hamiltonian["Hopping"])
        self.__MaxTauBin=self.__Map.MaxTauBin
        self.__Beta=self.__Map.Beta

        self.NearestNeighbor=[]
        for i in range(Lat.NSublat):
            self.NearestNeighbor.append([])
            for j in range(Lat.NSublat):
                self.NearestNeighbor[i].append([])

        self.NextNearestNeighbor=[]
        for i in range(Lat.NSublat):
            self.NextNearestNeighbor.append([])
            for j in range(Lat.NSublat):
                self.NextNearestNeighbor[i].append([])
        

    def Build(self):
        #self.BareG and self.BareW must be reinitialized at every time Build() is called
        self.BareG=weight.Weight("SmoothT", self.__Map, "TwoSpins", "AntiSymmetric", "R", "T")
        self.BareW=weight.Weight("DeltaT", self.__Map, "FourSpins", "Symmetric", "R", "T")
        LatName=self.Lat.Name
        if self.__Model=="Heisenberg":
            self.__Heisenberg(LatName)
        if self.__Model=="J1J2":
            self.__J1J2(LatName)
        elif self.__Model=="DiagCount":
            self.__DiagCount()
        return (self.BareG,self.BareW)

    def DecreaseField(self, Anneal):
        if abs(self.__DeltaField[0])>1e-5:
            for i in range(len(self.__DeltaField)):
                self.__DeltaField[i] += Anneal["Interval"][i]
                Anneal["DeltaField"][i] += Anneal["Interval"][i] 
        log.info(green("ExternalField decreased to: {0}".format(self.__DeltaField)))

    def RevertField(self, Anneal):
        for i in range(len(self.__DeltaField)):
            Anneal["Interval"][i]/=2.0
            self.__DeltaField[i] -= Anneal["Interval"][i]
            Anneal["DeltaField"][i] -= Anneal["Interval"][i]
        log.info(green("ExternalField reverted to: {0}".format(self.__DeltaField)))

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
                Mu=1j*np.pi/2.0/Beta+Pauli_Z[sp, sp]*(self.__DeltaField[sub]+self.__ExternalField[sub])
                self.BareG.Data[sp,sub,sp,sub,0,:]=np.exp(Mu*TauGrid)/(1.0+np.exp(Mu*Beta))

        Interaction=list(self.__Interaction)+[0,0,0,0,0]
        J1,J2=Interaction[0:2]
        #Bare W
        #Dimension: 2
        spin=self.__Map.Spin2Index(UP,UP)
        if LatName=="Checkerboard":
        #NSublat: 2
            Lx,Ly=self.__Map.L
            A,B=0,1
            self.NearestNeighbor[A][B]=[(0, 0),(0,Ly-1),(Lx-1,0),(Lx-1,Ly-1)]
            self.NearestNeighbor[B][A]=[(0, 0),(0,   1),(1,   0),(   1,   1)]

            self.NextNearestNeighbor[A][A]=[(0, 1),(1,   0),(0,Ly-1),(Lx-1,   0)]
            self.NextNearestNeighbor[B][B]=[(0, 1),(1,   0),(0,Ly-1),(Lx-1,   0)]

            #J1 interaction A-->B, B-->A
            for i in range(2):
                for j in range(2):
                    for e in self.NearestNeighbor[i][j]:
                        self.BareW.Data[spin,i,spin,j,self.__Map.CoordiIndex(e)] = J1/4;

            #J2 interaction A-->A, B-->B
            for i in range(2):
                for j in range(2):
                    for e in self.NextNearestNeighbor[i][j]:
                        self.BareW.Data[spin,i,spin,j,self.__Map.CoordiIndex(e)] = J2/4;
        elif LatName=="Square":
        #NSublat: 1
            Lx,Ly=self.__Map.L
            sub=0

            self.NearestNeighbor[0][0]=[(0,1),(1,0),(Lx-1,0),(0,Ly-1)]
            self.NextNearestNeighbor[0][0]=[(1,1),(Lx-1,1),(1,Ly-1),(Lx-1,Ly-1)]

            #J1 interaction on nearest neighbors
            for i in self.NearestNeighbor[0][0]:
                self.BareW.Data[spin,sub,spin,sub,self.__Map.CoordiIndex(i)] = J1/4.0;
            #J2 interaction on next nearest neighbors
            for i in self.NextNearestNeighbor[0][0]:
                self.BareW.Data[spin,sub,spin,sub,self.__Map.CoordiIndex(i)] = J2/4.0;

        elif LatName=="Honeycomb":
            #NSublat: 2
            Lx,Ly=self.__Map.L
            A,B=0,1
            self.NearestNeighbor[A][B]=[(0, 0), (Lx-1, 0), (Lx-1,Ly-1)]
            self.NearestNeighbor[B][A]=[(0, 0), (1,0), (1,1)]

            #J1 interaction A-->B, B-->A
            for i in range(2):
                for j in range(2):
                    for e in self.NearestNeighbor[i][j]:
                        self.BareW.Data[spin,i,spin,j,self.__Map.CoordiIndex(e)] = J1/4;
            ##J2 interaction A-->A, B-->B
        elif LatName=="Kagome":
            #NSublat: 2
            Lx,Ly=self.__Map.L
            A,B,C=0,1,2
            self.NearestNeighbor[A][B]=[(0, 0), (1, Ly-1)]
            self.NearestNeighbor[A][C]=[(0, 0), (0, Ly-1)]
            self.NearestNeighbor[B][A]=[(0, 0), (Lx-1, 1)]
            self.NearestNeighbor[B][C]=[(0, 0), (Lx-1, 0)]
            self.NearestNeighbor[C][A]=[(0, 0), (0, 1)]
            self.NearestNeighbor[C][B]=[(0, 0), (1, 0)]
            for i in range(3):
                for j in range(3):
                    for e in self.NearestNeighbor[i][j]:
                        self.BareW.Data[spin,i,spin,j,self.__Map.CoordiIndex(e)] = J1/4;

        elif LatName=="Cubic":
        #NSublat: 1
            Lx,Ly,Lz=self.__Map.L
            sub=0
            self.NearestNeighbor[0][0]=[(1,0,0),(0,1,0),(0,0,1),(Lx-1,0,0),(0,Ly-1,0),(0,0,Lz-1)]
            self.NextNearestNeighbor[0][0]=[]
            #J1 interaction on nearest neighbors
            for i in NearestNeighbor[0][0]:
                self.BareW.Data[spin,sub,spin,sub,self.__Map.CoordiIndex(i)] = J1/4.0;
            #J2 interaction on next nearest neighbors
            for i in self.NextNearestNeighbor[0][0]:
                self.BareW.Data[spin,sub,spin,sub,self.__Map.CoordiIndex(i)] = J2/4.0;
            pass

        elif LatName=="Pyrochlore":
        #NSublat: 4
            Lx,Ly,Lz=self.__Map.L
            A,B,C,D=0,1,2,3

            self.NearestNeighbor[A][B]=[(0,0,0),(0,0,Lz-1)]
            self.NearestNeighbor[A][C]=[(0,0,0),(Lx-1,0,0)]
            self.NearestNeighbor[A][D]=[(0,0,0),(0,Ly-1,0)]
            self.NearestNeighbor[B][A]=[(0,0,0),(0,0,1)]
            self.NearestNeighbor[B][C]=[(0,0,0),(Lx-1,0,1)]
            self.NearestNeighbor[B][D]=[(0,0,0),(0,Ly-1,1)]
            self.NearestNeighbor[C][A]=[(0,0,0),(1,0,0)]
            self.NearestNeighbor[C][B]=[(0,0,0),(1,0,Lz-1)]
            self.NearestNeighbor[C][D]=[(0,0,0),(1,Ly-1,0)]
            self.NearestNeighbor[D][A]=[(0,0,0),(0,1,0)]
            self.NearestNeighbor[D][B]=[(0,0,0),(0,1,Lz-1)]
            self.NearestNeighbor[D][C]=[(0,0,0),(Lx-1,1,0)]

            for i in range(4):
                for j in range(4):
                    for e in self.NearestNeighbor[i][j]:
                        self.BareW.Data[spin,i,spin,j,self.__Map.CoordiIndex(e)] = J1/4;

        elif LatName=="3DCheckerboard":
        #NSublat: 2
            Lx,Ly,Lz=self.__Map.L
            A,B=0,1

            self.NearestNeighbor[A][B]=[(0, 0, 0),(0,Ly-1,0),(Lx-1,0,0), \
                    (Lx-1,Ly-1,0),(Lx-1,Ly-1,1),(0,0,Lz-1)]
            self.NearestNeighbor[B][A]=[(0, 0, 0),(0,   1,0),(1,   0,0), \
                    (   1,   1,0),(0,0,1),(1,1,Lz-1)]

            self.NextNearestNeighbor[A][A]=[(0,1,0),(1,0,0),(0,Ly-1,0),(Lx-1,0,0),
                    (Lx-1,Ly-1,1),(0,0,Lz-1),(1,1,Lz-1),(0,0,1),
                    (0,Ly-1,1),(0,1,Lz-1),(1,0,Lz-1),(Lx-1,0,1)]

            self.NextNearestNeighbor[B][B]=[(0,1,0),(1,0,0),(0,Ly-1,0),(Lx-1,0,0),
                    (Lx-1,Ly-1,1),(1,1,Lz-1),(0,0,Lz-1),(0,0,1),
                    (1,0,Lz-1),(0,1,Lz-1),(Lx-1,0,1),(0,Ly-1,1)]

            #J1 interaction A-->B, B-->A
            for i in range(2):
                for j in range(2):
                    for e in self.NearestNeighbor[i][j]:
                        self.BareW.Data[spin,i,spin,j,self.__Map.CoordiIndex(e)] = J1/4;
            ##J2 interaction A-->A, B-->B
            for i in range(2):
                for j in range(2):
                    for e in self.NextNearestNeighbor[i][j]:
                        self.BareW.Data[spin,i,spin,j,self.__Map.CoordiIndex(e)] = J2/4;
        else:
            Assert(False, "Lattice {0} has not been implemented yet!".format(LatName))

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
        points, lines=self.Lat.GetSitesList()
        interaction=self.__GetBareWList()
        return {"Points":points, "Interaction":interaction, "Lines": lines}

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
        Spin=self.__Map.Spin2Index(UP,UP)
        offset=np.array(self.__Map.L)/2-1
        size=self.__Map.Vol*self.__Map.NSublat
        BareWList=[]
        Origin=[0 for e in self.__Map.L]
        for sub in self.__Map.GetAllSublatTuple():
            for coord in self.__Map.GetAllCoordi():
                weight=self.BareW.Data[Spin,sub[IN],Spin,sub[OUT],self.__Map.CoordiIndex(coord)]
                if weight*weight>1.0e-10:
                    n=self.__Map.LatIndex(coord, sub[OUT])
                    # vec is (real vector of in-site, real vector of out-site) 
                    vec=(self.Lat.GetRealVec(Origin,0, sub[IN], offset), \
                                    self.Lat.GetRealVec(coord,0, sub[OUT], offset))
                    coord=(self.__Map.LatIndex(Origin, sub[IN]), \
                                      self.__Map.LatIndex(coord, sub[OUT]))
                    BareWList.append([vec, coord, sub[IN]])
        return BareWList

