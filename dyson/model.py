import weight
import lattice as lat
from weight import UP,DOWN,IN,OUT
import numpy as np
import math
import traceback
from logger import *

class BareFactory:
    def __init__(self, _map, Lat, Hamiltonian, Anneal):
        self.Lat=Lat
        self.__Map=_map
        self.__Model=Hamiltonian["Name"]
        self.__Interaction=np.array(Hamiltonian["Interaction"])
        self.__ExternalField=np.array(Hamiltonian["ExternalField"])
        self.__DeltaField=np.array(Anneal["DeltaField"])
        self.__MaxTauBin=self.__Map.MaxTauBin
        self.__Beta=self.__Map.Beta
        if "Hopping" in Hamiltonian:
            self.__Hopping=np.array(Hamiltonian["Hopping"])
        if "Phase" in Hamiltonian:
            self.__Phase=np.array(Hamiltonian["Phase"])
        if "HubbardInteraction" in Hamiltonian:
            self.__HubbardInteraction=np.array(Hamiltonian["HubbardInteraction"])
        if "ShortRangeInteraction" in Hamiltonian:
            self.__ShortRangeInteraction=np.array(Hamiltonian["ShortRangeInteraction"])
        if "CoulombInteraction" in Hamiltonian:
            self.__ShortRangeInteraction=np.array(Hamiltonian["CoulombInteraction"])
        if "ChemicalPotential" in Hamiltonian:
            self.__Mu=np.array(Hamiltonian["ChemicalPotential"])
        if "Description" in Hamiltonian:
            self.__Description=Hamiltonian["Description"]
        else:
            self.__Description=None

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

        self.NextNextNearestNeighbor=[]
        for i in range(Lat.NSublat):
            self.NextNextNearestNeighbor.append([])
            for j in range(Lat.NSublat):
                self.NextNextNearestNeighbor[i].append([])
        

    def Build(self):
        #self.BareG and self.BareW must be reinitialized at every time Build() is called
        self.BareG=weight.Weight("SmoothT", self.__Map, "TwoSpins", "AntiSymmetric", "R", "T")
        self.BareW=weight.Weight("DeltaT", self.__Map, "FourSpins", "Symmetric", "R", "T")
        LatName=self.Lat.Name
        # getattr(self, self.__Model)(LatName)
        try:
            getattr(self, self.__Model)(LatName)
        except:
            log.error(blue("Model construction fails {0}".format(traceback.format_exc())))
            # Assert(False, "Model {0} has not been implemented!".format(self.__Model))

        return (self.BareG,self.BareW)

    def DecreaseField(self, Anneal):
        #TODO: what if DeltaField/Interval is not an integer?!
        flag=False
        if abs(self.__DeltaField[0])>1e-5:
            for i in range(len(self.__DeltaField)):
                self.__DeltaField[i] += Anneal["Interval"][i]
                Anneal["DeltaField"][i] += Anneal["Interval"][i] 
            flag=True
        log.info(green("ExternalField decreased to: {0}".format(self.__DeltaField)))
        return flag

    def RevertField(self, Anneal):
        for i in range(len(self.__DeltaField)):
            Anneal["Interval"][i]/=2.0
            self.__DeltaField[i] -= Anneal["Interval"][i]
            Anneal["DeltaField"][i] -= Anneal["Interval"][i]
        log.info(green("ExternalField reverted to: {0}".format(self.__DeltaField)))

    #model defintion
    def DiagCount(self, LatName):
        raise NotImplementedError
    def Hubbard(self, LatName):
        self.__BuildBareG(LatName)
        self.__BuildHubbardInteraction()
    def Haldane(self, LatName):
        Assert(LatName=="Honeycomb", "Haldane model is defined on the Honeycomb lattice!")
        self.__BuildBareG(LatName="Honeycomb")
        self.__BuildShortRangeInteraction(LatName="Honeycomb", InteractionTypeList=["nn"])
    def Heisenberg(self, LatName):
        Assert(len(self.__Interaction)>=1, "Heisenberg model only has one coupling!")
        self.__SpinModel(LatName)
    def J1J2(self, LatName):
        Assert(len(self.__Interaction)>=2, "J1J2 model takes at least two couplings!")
        self.__SpinModel(LatName)
    def Kitaev(self, LatName):
        Assert(LatName=="Honeycomb", "Kitaev model takes three couplings!")
        self.__SpinModel("Kitaev")

    def __BuildBareG(self, LatName):
        Beta=self.__Map.Beta
        self.BareG.Data=np.zeros(self.BareG.Shape, dtype=complex)
        TauGrid=np.array([self.__Map.IndexToTau(t) for t in range(self.__MaxTauBin)])
        self.__KineticTerm=weight.Weight("DeltaT", self.__Map, "TwoSpins", "Symmetric", "R", "T")
        self.__LocalTerm=weight.Weight("DeltaT", self.__Map, "TwoSpins", "Symmetric", "R", "T")

        TauGrid=np.array([self.__Map.IndexToTau(t) for t in range(self.__MaxTauBin)])

        Hopping=list(self.__Hopping)+[0,0,0,0,0]
        Mu=list(self.__Mu)+[0,0,0,0,0]
        ExternalField=list(self.__ExternalField)+[0,0,0,0,0]
        Phase=list(self.__Phase)

        Pauli=self.__Map.Pauli()
        SpinOrbit=[]
        I=np.array([[1.0,0.0],[0.0,1.0]])
        for i in range(len(Hopping)):
            Phase=np.array(Phase)
            PhaseMod=np.sqrt(sum(Phase**2))
            if i<len(Phase) and PhaseMod>1e-8:
                # SpinOrbit.append(np.exp(-1j*(Phase[0]*Pauli[0]+Phase[1]*Pauli[1]+Phase[2]*Pauli[2])))
                PhaseMatrix=Phase[0]*Pauli[0]+Phase[1]*Pauli[1]+Phase[2]*Pauli[2]
                SpinOrbit.append(np.cos(PhaseMod)*I+1j*(PhaseMatrix/PhaseMod)*np.sin(PhaseMod))
            else:
                SpinOrbit.append(I)

        if LatName=="Square":
            #NSublat: 1
            Lx,Ly=self.__Map.L
            sub=0

            self.NearestNeighbor[0][0]=[(0,1),(1,0),(Lx-1,0),(0,Ly-1)]
            self.NextNearestNeighbor[0][0]=[(1,1),(Lx-1,1),(1,Ly-1),(Lx-1,Ly-1)]

            for e in self.NearestNeighbor[sub][sub]:
                self.__KineticTerm.Data[:,sub,:,sub,self.__Map.CoordiIndex(e)]+= -self.__Hopping[0]*SpinOrbit[0];
            for e in self.NextNearestNeighbor[sub][sub]:
                self.__KineticTerm.Data[:,sub,:,sub,self.__Map.CoordiIndex(e)]+= -self.__Hopping[1]*SpinOrbit[1];

        elif LatName=="Honeycomb":
            #NSublat: 2
            Lx,Ly=self.__Map.L
            A,B=0,1
            self.NearestNeighbor[A][B]=[(0, 0), (Lx-1, 0), (Lx-1,Ly-1)]
            self.NearestNeighbor[B][A]=[(0, 0), (1,0), (1,1)]

            self.NextNearestNeighbor[A][A]=[(1, 0), (Lx-1, 0), (0, 1), (0, Ly-1), (1, 1), (Lx-1,Ly-1)]
            self.NextNearestNeighbor[B][B]=[(1, 0), (Lx-1, 0), (0, 1), (0, Ly-1), (1, 1), (Lx-1,Ly-1)]

            #hopping A-->B, B-->A
            for i in range(2):
                for j in range(2):
                    for e in self.NearestNeighbor[i][j]:
                        self.__KineticTerm.Data[:,i,:,j,self.__Map.CoordiIndex(e)]+= -self.__Hopping[0]*SpinOrbit[0];
                    for e in self.NextNearestNeighbor[i][j]:
                        self.__KineticTerm.Data[:,i,:,j,self.__Map.CoordiIndex(e)]+= -self.__Hopping[1]*SpinOrbit[1];

        else:
            Assert(False, "Lattice {0} has not been implemented yet!".format(LatName))

        origin=self.__Map.CoordiIndex((0,0))
        for sub in range(self.__Map.NSublat):
            self.__LocalTerm.Data[UP,sub,UP,sub,origin]+= -Mu[0];
            self.__LocalTerm.Data[DOWN,sub,DOWN,sub,origin]+= -Mu[0];
            self.__LocalTerm.Data[UP,sub,UP,sub,origin]+= -ExternalField[sub];
            self.__LocalTerm.Data[DOWN,sub,DOWN,sub,origin]+= ExternalField[sub];

        # print "Kinectic in R:", self.__KineticTerm.Data[0,0,0,0,:]

        self.__KineticTerm.FFT("K")
        self.__LocalTerm.FFT("K")

        # print "Kinectic in K:", self.__KineticTerm.Data[0,0,0,0,:]

        self.__KineticTerm.Data+=self.__LocalTerm.Data
        Sp, Sub = self.__KineticTerm.NSpin, self.__KineticTerm.NSublat
        self.__KineticTerm.Data = self.__KineticTerm.Data.reshape([Sp*Sub, Sp*Sub, self.__Map.Vol])

        self.BareG.FFT("K","T")

        for k in range(self.__Map.Vol): 
            Ek,Uk=np.linalg.eig(self.__KineticTerm.Data[:,:,k])
            # print Ek, self.__KineticTerm.Data[:,:,k]
            Ukdag=Uk.conj().T
            for t in range(self.__MaxTauBin):
                Gk=np.zeros([Sp*Sub, Sp*Sub], dtype="complex")
                for i in range(Sp*Sub):
                    # Gk[i,i]=-np.exp(-Ek[i]*TauGrid[t])*(1.0-1.0/(1.0+np.exp(Beta*Ek[i])))
                    # Gk[i,i]=-np.exp(-Ek[i]*TauGrid[t])/(1.0+np.exp(-Beta*Ek[i]))
                    Gk[i,i]=self.__FermiKernel(Ek[i], TauGrid[t], Beta)
                Gk=np.dot(np.dot(Uk,Gk),Ukdag)
                self.BareG.Data[:,:,:,:,k,t]=Gk.reshape([Sp, Sub, Sp, Sub])

        # print "Gt:", self.BareG.Data[0,0,0,0,0,:]

    def __FermiKernel(self, w, t, beta):
        x=beta*w/2
        y=2*t/beta-1
        if x>100:
            return np.exp(-x*(y+1.))
        if x<-100:
            return np.exp(x*(1.0-y))
        return np.exp(-x*y)/(2*np.cosh(x))

    def __BuildHubbardInteraction(self):
        HubbardU=self.__HubbardInteraction
        for sub in range(self.__Map.NSublat):
            self.BareW.Data[0,sub,1,sub,:]+= HubbardU;
            self.BareW.Data[1,sub,0,sub,:]+= HubbardU;

    def __BuildShortRangeInteraction(self, LatName, InteractionTypeList):
        """
        InteractionTypeList: a list of types of interaction. For example, ["SS","SxSx","SzSz"] mean
        the Hamiltonian has three different interaction terms: Heisenberg term, SxSx term and SzSz term.
        """
        Beta=self.__Map.Beta
        self.BareW.Data=np.zeros(self.BareW.Shape, dtype=complex)
        Sx=0.5*np.array([0,1,1,0])
        Sy=0.5*np.array([0,-1j,1j,0])
        Sz=0.5*np.array([1,0,0,-1])
        I=np.array([1,0,0,1])
        SS=np.outer(Sx,Sx)+np.outer(Sy,Sy)+np.outer(Sz,Sz)
        SSxy=np.outer(Sx,Sx)+np.outer(Sy,Sy)
        SxSx=np.outer(Sx,Sx)
        SySy=np.outer(Sy,Sy)
        SzSz=np.outer(Sz,Sz)
        II=np.outer(I,I)

        ShortRangeInteraction=list(self.__ShortRangeInteraction)+[0,0,0,0,0]

        TypeList=[]
        for t in InteractionTypeList:
            if t=="SS":
                TypeList.append(SS)
            elif t=="SSxy":
                TypeList.append(SSxy)
            elif t=="SxSx":
                TypeList.append(SxSx)
            elif t=="SySy":
                TypeList.append(SySy)
            elif t=="SzSz":
                TypeList.append(SzSz)
            elif t=="nn":
                TypeList.append(II)
        for i in range(len(TypeList), len(ShortRangeInteraction)):
            TypeList.append(0.0)

        if LatName=="Square":
        #NSublat: 1
            Lx,Ly=self.__Map.L
            sub=0

            self.NearestNeighbor[0][0]=[(0,1),(1,0),(Lx-1,0),(0,Ly-1)]
            self.NextNearestNeighbor[0][0]=[(1,1),(Lx-1,1),(1,Ly-1),(Lx-1,Ly-1)]

            #J1 interaction on nearest neighbors
            for i in self.NearestNeighbor[0][0]:
                self.BareW.Data[:,sub,:,sub,self.__Map.CoordiIndex(i)]+= ShortRangeInteraction[0]*TypeList[0];
            #J2 interaction on next nearest neighbors
            for i in self.NextNearestNeighbor[0][0]:
                self.BareW.Data[:,sub,:,sub,self.__Map.CoordiIndex(i)]+= ShortRangeInteraction[1]*TypeList[1];

        elif LatName=="Honeycomb":
            #NSublat: 2
            Lx,Ly=self.__Map.L
            A,B=0,1
            self.NearestNeighbor[A][B]=[(0, 0), (Lx-1, 0), (Lx-1,Ly-1)]
            self.NearestNeighbor[B][A]=[(0, 0), (1,0), (1,1)]

            self.NextNearestNeighbor[A][A]=[(1, 0), (Lx-1, 0), (0, 1), (0, Ly-1), (1, 1), (Lx-1,Ly-1)]
            self.NextNearestNeighbor[B][B]=[(1, 0), (Lx-1, 0), (0, 1), (0, Ly-1), (1, 1), (Lx-1,Ly-1)]

            #J1 interaction A-->B, B-->A
            for i in range(2):
                for j in range(2):
                    for e in self.NearestNeighbor[i][j]:
                        self.BareW.Data[:,i,:,j,self.__Map.CoordiIndex(e)]+= ShortRangeInteraction[0]*TypeList[0];
                        # self.BareW.Data[:,i,:,j,self.__Map.CoordiIndex(e)]+= J1*SSxy;

                    ##J2 interaction A-->A, B-->B
                    for e in self.NextNearestNeighbor[i][j]:
                        self.BareW.Data[:,i,:,j,self.__Map.CoordiIndex(e)]+= ShortRangeInteraction[1]*TypeList[1];
                        # self.BareW.Data[:,i,:,j,self.__Map.CoordiIndex(e)]+= J2*SSxy;
        else:
            Assert(False, "Lattice {0} has not been implemented yet!".format(LatName))

    def __SpinModel(self, LatName):
        Beta=self.__Map.Beta
        self.BareG.Data=np.zeros(self.BareG.Shape, dtype=complex)
        self.BareW.Data=np.zeros(self.BareW.Shape, dtype=complex)
        Sx=0.5*np.array([0,1,1,0])
        Sy=0.5*np.array([0,-1j,1j,0])
        Sz=0.5*np.array([1,0,0,-1])
        I=np.array([1,0,0,1])
        SS=np.outer(Sx,Sx)+np.outer(Sy,Sy)+np.outer(Sz,Sz)
        SSxy=np.outer(Sx,Sx)+np.outer(Sy,Sy)
        SxSx=np.outer(Sx,Sx)
        SySy=np.outer(Sy,Sy)
        SzSz=np.outer(Sz,Sz)

        if self.__Description is not None and "ImW" in self.__Description:
            #use imaginary W instead of imaginary chemical potential
            II=np.outer(I,I)
            #for i in range(self.__Map.NSublat):
                #self.BareW.Data[:,i,:,i,0]+=1j*np.pi/4.0*II/Beta
            #self.__Mu=1j*np.pi/4.0/Beta
            self.__Mu=np.pi/4.0/Beta
        else:
            self.__Mu=1j*np.pi/2.0/Beta
        #Bare G
        self.__Hopping=np.array([0.0])
        log.info("set Mu={0}, Hopping={1}, and SmoothT Bare G".format(self.__Mu, self.__Hopping))
        Assert(len(self.__ExternalField)>=self.__Map.NSublat, 
                "expect at least {0} externalfield components!".format(self.__Map.NSublat))

        TauGrid=np.array([self.__Map.IndexToTau(t) for t in range(self.__MaxTauBin)])
        Pauli_Z=self.__Map.Pauli()[2]
        for sub in range(self.__Map.NSublat):
            for sp in range(2):
                Mu=self.__Mu+0.5*Pauli_Z[sp, sp]*(self.__DeltaField[sub]+self.__ExternalField[sub])
                self.BareG.Data[sp,sub,sp,sub,0,:]=np.exp(Mu*TauGrid)/(1.0+np.exp(Mu*Beta))

        Interaction=list(self.__Interaction)+[0,0,0,0,0]
        J1,J2,J3=Interaction[0:3]
        J_perturbation=None
        if len(Interaction)>3:
            J_perturbation=Interaction[3:]
        #Bare W
        #Dimension: 2
        spin=self.__Map.Spin2Index(UP,UP)

        if LatName=="Chain":
        #NSublat: 1
            Lx=self.__Map.L[0]  #for 1D system, L=[Lx]
            sub=0

            self.NearestNeighbor[0][0]=[(1,),(Lx-1,)]
            if Lx>2:
                self.NextNearestNeighbor[0][0]=[(2,),(Lx-2,)]
            elif Lx==2:
                self.NextNearestNeighbor[0][0]=[(0,),(0,)]
            else:
                raise ValueError("System size must be larger than 1")

            #J1 interaction on nearest neighbors
            for i in self.NearestNeighbor[0][0]:
                self.BareW.Data[:,sub,:,sub,self.__Map.CoordiIndex(i)]+= J1*SS;
            #J2 interaction on next nearest neighbors
            for i in self.NextNearestNeighbor[0][0]:
                self.BareW.Data[:,sub,:,sub,self.__Map.CoordiIndex(i)]+= J2*SS;

        elif LatName=="Checkerboard":
        #NSublat: 2
            Lx,Ly=self.__Map.L
            A,B=0,1
            self.NearestNeighbor[A][B]=[(0, 0),(0,Ly-1),(Lx-1,0),(Lx-1,Ly-1)]
            self.NearestNeighbor[B][A]=[(0, 0),(0,   1),(1,   0),(   1,   1)]

            self.NextNearestNeighbor[A][A]=[(0, 1),(1,   0),(0,Ly-1),(Lx-1,   0)]
            self.NextNearestNeighbor[B][B]=[(0, 1),(1,   0),(0,Ly-1),(Lx-1,   0)]

            for i in range(2):
                for j in range(2):
                    #J1 interaction A-->B, B-->A
                    for e in self.NearestNeighbor[i][j]:
                        self.BareW.Data[:,i,:,j,self.__Map.CoordiIndex(e)]+= J1*SS;
                        #self.BareW.Data[:,i,:,j,self.__Map.CoordiIndex(e)]+= J1*SzSz;
                    #J2 interaction A-->A, B-->B
                    for e in self.NextNearestNeighbor[i][j]:
                        self.BareW.Data[:,i,:,j,self.__Map.CoordiIndex(e)]+= J2*SS;

        elif LatName=="ValenceBond":
            Lx,Ly=self.__Map.L
            A,B=0,1
            
            self.NearestNeighbor[A][B]=[(0, 0)]
            self.NearestNeighbor[B][A]=[(0, 0)]
            self.NearestNeighbor[A][A]=[(Lx-1, 0),(1,0)]
            self.NearestNeighbor[B][B]=[(Lx-1, 0),(1,0)]

            self.NextNearestNeighbor[A][B]=[(0,   Ly-1)]
            self.NextNearestNeighbor[B][A]=[(0,   1)]

            # self.NearestNeighbor[A][B]=[(0, 0),(0,Ly-1)]
            # self.NearestNeighbor[B][A]=[(0, 0),(0,   1)]
            # self.NearestNeighbor[A][A]=[(Lx-1, 0),(1,0)]
            # self.NearestNeighbor[B][B]=[(Lx-1, 0),(1,0)]

            # self.NextNearestNeighbor[A][B]=[(1, 0),(1,   Ly-1),(Lx-1,0),(Lx-1,   Ly-1)]
            # self.NextNearestNeighbor[B][A]=[(1, 1),(1,   0),(Lx-1,1),(Lx-1,   0)]

            for i in range(2):
                for j in range(2):
                    #J1 interaction A-->B, B-->A
                    for e in self.NearestNeighbor[i][j]:
                        self.BareW.Data[:,i,:,j,self.__Map.CoordiIndex(e)]+= J1*SS;
                    #J2 interaction A-->A, B-->B
                    for e in self.NextNearestNeighbor[i][j]:
                        self.BareW.Data[:,i,:,j,self.__Map.CoordiIndex(e)]+= J2*SS;
            self.BareW.Data[:,A,:,B,self.__Map.CoordiIndex((0,0))]+= J_perturbation[0]*SS;
            self.BareW.Data[:,B,:,A,self.__Map.CoordiIndex((0,0))]+= J_perturbation[0]*SS;

        elif LatName=="Square":
        #NSublat: 1
            Lx,Ly=self.__Map.L
            sub=0

            self.NearestNeighbor[0][0]=[(0,1),(1,0),(Lx-1,0),(0,Ly-1)]
            self.NextNearestNeighbor[0][0]=[(1,1),(Lx-1,1),(1,Ly-1),(Lx-1,Ly-1)]

            #J1 interaction on nearest neighbors
            for i in self.NearestNeighbor[0][0]:
                self.BareW.Data[:,sub,:,sub,self.__Map.CoordiIndex(i)]+= J1*SS;
            #J2 interaction on next nearest neighbors
            for i in self.NextNearestNeighbor[0][0]:
                self.BareW.Data[:,sub,:,sub,self.__Map.CoordiIndex(i)]+= J2*SS;

        elif LatName=="Triangular":
        #NSublat: 1
            Lx,Ly=self.__Map.L
            sub=0

            self.NearestNeighbor[0][0]=[(0,1),(1,0),(1,Ly-1),(0,Ly-1), (Lx-1,0), (Lx-1, 1)]
            self.NextNearestNeighbor[0][0]=[(1,1),(2,Ly-1),(1,Ly-2),(Lx-1,Ly-1),(Lx-2,1), (Lx-1,2)]

            #J1 interaction on nearest neighbors
            for i in self.NearestNeighbor[0][0]:
                self.BareW.Data[:,sub,:,sub,self.__Map.CoordiIndex(i)]+= J1*SS;
            #J2 interaction on next nearest neighbors
            for i in self.NextNearestNeighbor[0][0]:
                self.BareW.Data[:,sub,:,sub,self.__Map.CoordiIndex(i)]+= J2*SS;

        elif LatName=="Assymetric_Triangular":
        #NSublat: 3
            Lx,Ly=self.__Map.L

            A,B,C=0,1,2

            self.NearestNeighbor[A][B]=[(0, 0), (0, Ly-1), (Lx-1, 0)]
            self.NearestNeighbor[A][C]=[(0, 0), (Lx-1, 1), (Lx-1, 0)]
            self.NearestNeighbor[B][A]=[(0, 0), (0, 1), (1, 0)]
            self.NearestNeighbor[B][C]=[(0, 0), (0, 1), (Lx-1, 1)]
            self.NearestNeighbor[C][A]=[(0, 0), (1, 0), (1, Ly-1)]
            self.NearestNeighbor[C][B]=[(0, 0), (1, Ly-1), (0, Ly-1)]

            for i in range(3):
                for j in range(3):
                    for e in self.NearestNeighbor[i][j]:
                        self.BareW.Data[:,i,:,j,self.__Map.CoordiIndex(e)]+= J1*SS;

        elif LatName=="Honeycomb":
            #NSublat: 2
            Lx,Ly=self.__Map.L
            A,B=0,1
            self.NearestNeighbor[A][B]=[(0, 0), (Lx-1, 0), (Lx-1,Ly-1)]
            self.NearestNeighbor[B][A]=[(0, 0), (1,0), (1,1)]

            self.NextNearestNeighbor[A][A]=[(1, 0), (Lx-1, 0), (0, 1), (0, Ly-1), (1, 1), (Lx-1,Ly-1)]
            self.NextNearestNeighbor[B][B]=[(1, 0), (Lx-1, 0), (0, 1), (0, Ly-1), (1, 1), (Lx-1,Ly-1)]

            #J1 interaction A-->B, B-->A
            for i in range(2):
                for j in range(2):
                    for e in self.NearestNeighbor[i][j]:
                        self.BareW.Data[:,i,:,j,self.__Map.CoordiIndex(e)]+= J1*SS;
                        # self.BareW.Data[:,i,:,j,self.__Map.CoordiIndex(e)]+= J1*SSxy;

                    ##J2 interaction A-->A, B-->B
                    for e in self.NextNearestNeighbor[i][j]:
                        self.BareW.Data[:,i,:,j,self.__Map.CoordiIndex(e)]+= J2*SS;
                        # self.BareW.Data[:,i,:,j,self.__Map.CoordiIndex(e)]+= J2*SSxy;

        elif LatName=="Kitaev":
            #NSublat: 2
            Lx,Ly=self.__Map.L
            A,B=0,1

            self.BareW.Data[:,A,:,B,self.__Map.CoordiIndex((0,0))]+= J1*SxSx;
            self.BareW.Data[:,B,:,A,self.__Map.CoordiIndex((0,0))]+= J1*SxSx;

            self.BareW.Data[:,A,:,B,self.__Map.CoordiIndex((Lx-1,0))]+= J2*SySy;
            self.BareW.Data[:,B,:,A,self.__Map.CoordiIndex((1,0))]+= J2*SySy;

            self.BareW.Data[:,A,:,B,self.__Map.CoordiIndex((Lx-1,Ly-1))]+= J3*SzSz;
            self.BareW.Data[:,B,:,A,self.__Map.CoordiIndex((1,1))]+= J3*SzSz;

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

            self.NextNearestNeighbor[A][B]=[(1, 0), (0, Ly-1)]
            self.NextNearestNeighbor[A][C]=[(1, Ly-1), (Lx-1, 0)]
            self.NextNearestNeighbor[B][A]=[(0, 1), (Lx-1, 0)]
            self.NextNearestNeighbor[B][C]=[(Lx-1, 1), (0, Ly-1)]
            self.NextNearestNeighbor[C][A]=[(1, 0), (Lx-1, 1)]
            self.NextNearestNeighbor[C][B]=[(0, 1), (1, Ly-1)]

            self.NextNextNearestNeighbor[A][A]=[(1, 0), (Lx-1,0)]
            self.NextNextNearestNeighbor[B][B]=[(0, 1), (0,Ly-1)]
            self.NextNextNearestNeighbor[C][C]=[(Lx-1, 1), (1,Ly-1)]

            for i in range(3):
                for j in range(3):
                    #J1 interaction
                    for e in self.NearestNeighbor[i][j]:
                        self.BareW.Data[:,i,:,j,self.__Map.CoordiIndex(e)]+= J1*SS;

                    ##J2 interaction
                    for e in self.NextNearestNeighbor[i][j]:
                        self.BareW.Data[:,i,:,j,self.__Map.CoordiIndex(e)]+= J2*SS;
                        # self.BareW.Data[:,i,:,j,self.__Map.CoordiIndex(e)]+= J2*SSxy;

                    ##J3 interaction
                    for e in self.NextNextNearestNeighbor[i][j]:
                        self.BareW.Data[:,i,:,j,self.__Map.CoordiIndex(e)]+= J3*SS;
                        # self.BareW.Data[:,i,:,j,self.__Map.CoordiIndex(e)]+= J2*SSxy;

        elif LatName=="Cubic":
        #NSublat: 1
            Lx,Ly,Lz=self.__Map.L
            sub=0
            self.NearestNeighbor[0][0]=[(1,0,0),(0,1,0),(0,0,1),(Lx-1,0,0),(0,Ly-1,0),(0,0,Lz-1)]
            self.NextNearestNeighbor[0][0]=[]
            #J1 interaction on nearest neighbors
            for i in self.NearestNeighbor[0][0]:
                self.BareW.Data[:,sub,:,sub,self.__Map.CoordiIndex(i)]+= J1*SS;
            #J2 interaction on next nearest neighbors
            for i in self.NextNearestNeighbor[0][0]:
                self.BareW.Data[:,sub,:,sub,self.__Map.CoordiIndex(i)]+= J2*SS;

        elif LatName=="Pyrochlore":
        #NSublat: 4
            Lx,Ly,Lz=self.__Map.L
            A,B,C,D=0,1,2,3

            self.NearestNeighbor[A][B]=[(0,0,0),(Lx-1,0,0)]
            self.NearestNeighbor[A][C]=[(0,0,0),(0,Ly-1,0)]
            self.NearestNeighbor[A][D]=[(0,0,0),(0,0,Lz-1)]
            self.NearestNeighbor[B][A]=[(0,0,0),(1,0,0)]
            self.NearestNeighbor[B][C]=[(0,0,0),(1,Ly-1,0)]
            self.NearestNeighbor[B][D]=[(0,0,0),(1,0,Lz-1)]
            self.NearestNeighbor[C][A]=[(0,0,0),(0,1,0)]
            self.NearestNeighbor[C][B]=[(0,0,0),(Lx-1,1,0)]
            self.NearestNeighbor[C][D]=[(0,0,0),(0,1,Lz-1)]
            self.NearestNeighbor[D][A]=[(0,0,0),(0,0,1)]
            self.NearestNeighbor[D][B]=[(0,0,0),(Lx-1,0,1)]
            self.NearestNeighbor[D][C]=[(0,0,0),(0,Ly-1,1)]

            self.NextNearestNeighbor[A][B]=[(0,Ly-1,0),(0,0,Lz-1),(Lx-1,1,0),(Lx-1,0,1)]
            self.NextNearestNeighbor[A][C]=[(Lx-1,0,0),(0,0,Lz-1),(1,Ly-1,0),(0,Ly-1,1)]
            self.NextNearestNeighbor[A][D]=[(Lx-1,0,0),(0,Ly-1,0),(1,0,Lz-1),(0,1,Lz-1)]
            self.NextNearestNeighbor[B][A]=[(1,Ly-1,0),(1,0,Lz-1),(0,0,1),(0,1,0)]
            self.NextNearestNeighbor[B][C]=[(1,0,0),(1,0,Lz-1),(0,Ly-1,0),(0,Ly-1,1)]
            self.NextNearestNeighbor[B][D]=[(1,Ly-1,0),(1,0,0),(0,0,Lz-1),(0,1,Lz-1)]
            self.NextNearestNeighbor[C][A]=[(Lx-1,1,0),(0,1,Lz-1),(1,0,0),(0,0,1)]
            self.NextNearestNeighbor[C][B]=[(0,1,0),(0,1,Lz-1),(Lx-1,0,0),(Lx-1,0,1)]
            self.NextNearestNeighbor[C][D]=[(0,1,0),(Lx-1,1,0),(0,0,Lz-1),(1,0,Lz-1)]
            self.NextNearestNeighbor[D][A]=[(Lx-1,0,1),(0,Ly-1,1),(1,0,0),(0,1,0)]
            self.NextNearestNeighbor[D][B]=[(0,0,1),(0,Ly-1,1),(Lx-1,0,0),(Lx-1,1,0)]
            self.NextNearestNeighbor[D][C]=[(0,0,1),(Lx-1,0,1),(0,Ly-1,0),(1,Ly-1,0),]

            ########J3, only connected####################
            self.NextNextNearestNeighbor[A][A]=[(1,0,0),(0,1,0),(0,0,1),(Lx-1,0,0),(0,Ly-1,0),(0,0,Lz-1)]
            self.NextNextNearestNeighbor[B][B]=[(Lx-1,0,0),(Lx-1,1,0),(Lx-1,0,1),(1,0,0),(1,Ly-1,0),(1,0,Lz-1)]
            self.NextNextNearestNeighbor[C][C]=[(0,Ly-1,0),(1,Ly-1,0),(0,Ly-1,1),(0,1,0),(Lx-1,1,0),(0,1,Lz-1)]
            self.NextNextNearestNeighbor[D][D]=[(0,0,Lz-1),(1,0,Lz-1),(0,1,Lz-1),(0,0,1),(Lx-1,0,1),(0,Ly-1,1)]

            for i in range(4):
                for j in range(4):
                    for e in self.NearestNeighbor[i][j]:
                        #self.BareW.Data[:,i,:,j,self.__Map.CoordiIndex(e)]+= J1*SzSz;
                        self.BareW.Data[:,i,:,j,self.__Map.CoordiIndex(e)]+= J1*SS;
                    for e in self.NextNearestNeighbor[i][j]:
                        #self.BareW.Data[:,i,:,j,self.__Map.CoordiIndex(e)]+= J1*SzSz;
                        self.BareW.Data[:,i,:,j,self.__Map.CoordiIndex(e)]+= J2*SS;
                    for e in self.NextNextNearestNeighbor[i][j]:
                        #self.BareW.Data[:,i,:,j,self.__Map.CoordiIndex(e)]+= J1*SzSz;
                        self.BareW.Data[:,i,:,j,self.__Map.CoordiIndex(e)]+= J3*SS;
            #S111S111=np.outer(Sx+Sy+Sz,Sx+Sy+Sz)
            #for i in range(4):
                #self.BareW.Data[:,i,:,i,0]+ = -10.0*J1*S111S111;

        elif LatName=="DecoPyrochlore":
        #NSublat: 6
            Lx,Ly,Lz=self.__Map.L
            A,B,C,D,E,F=0,1,2,3,4,5

            # self.NearestNeighbor[A][B]=[(0,0,0),(Lx-1,0,0)]
            # self.NearestNeighbor[A][C]=[(0,0,0),(0,Ly-1,0)]
            # self.NearestNeighbor[A][D]=[(0,0,0),(0,0,Lz-1)]
            # self.NearestNeighbor[B][A]=[(0,0,0),(1,0,0)]
            # self.NearestNeighbor[B][C]=[(0,0,0),(1,Ly-1,0)]
            # self.NearestNeighbor[B][D]=[(0,0,0),(1,0,Lz-1)]
            # self.NearestNeighbor[C][A]=[(0,0,0),(0,1,0)]
            # self.NearestNeighbor[C][B]=[(0,0,0),(Lx-1,1,0)]
            # self.NearestNeighbor[C][D]=[(0,0,0),(0,1,Lz-1)]
            # self.NearestNeighbor[D][A]=[(0,0,0),(0,0,1)]
            # self.NearestNeighbor[D][B]=[(0,0,0),(Lx-1,0,1)]
            # self.NearestNeighbor[D][C]=[(0,0,0),(0,Ly-1,1)]

            self.NearestNeighbor[A][E]=[(0,0,0)]
            self.NearestNeighbor[B][E]=[(0,0,0)]
            self.NearestNeighbor[C][E]=[(0,0,0)]
            self.NearestNeighbor[D][E]=[(0,0,0)]
            self.NearestNeighbor[D][F]=[(0,0,0)]
            self.NearestNeighbor[E][A]=[(0,0,0)]
            self.NearestNeighbor[E][B]=[(0,0,0)]
            self.NearestNeighbor[E][C]=[(0,0,0)]
            self.NearestNeighbor[E][D]=[(0,0,0)]
            self.NearestNeighbor[F][D]=[(0,0,0)]
            self.NearestNeighbor[A][F]=[(Lx-1,0,0)]
            self.NearestNeighbor[B][F]=[(0,Ly-1,0)]
            self.NearestNeighbor[C][F]=[(0,0,Lz-1)]
            self.NearestNeighbor[F][A]=[(1,0,0)]
            self.NearestNeighbor[F][B]=[(0,1,0)]
            self.NearestNeighbor[F][C]=[(0,0,1)]

            self.NextNearestNeighbor[A][B]=[(0,0,0),(Lx-1,1,0)]
            self.NextNearestNeighbor[A][C]=[(0,0,0),(Lx-1,0,1)]
            self.NextNearestNeighbor[A][D]=[(0,0,0),(Lx-1,0,0)]
            self.NextNearestNeighbor[B][A]=[(0,0,0),(1,Ly-1,0)]
            self.NextNearestNeighbor[B][C]=[(0,0,0),(0,Ly-1,1)]
            self.NextNearestNeighbor[B][D]=[(0,0,0),(0,Ly-1,0)]
            self.NextNearestNeighbor[C][A]=[(0,0,0),(1,0,Lz-1)]
            self.NextNearestNeighbor[C][B]=[(0,0,0),(0,1,Lz-1)]
            self.NextNearestNeighbor[C][D]=[(0,0,0),(0,0,Lz-1)]
            self.NextNearestNeighbor[D][A]=[(0,0,0),(1,0,0)]
            self.NextNearestNeighbor[D][B]=[(0,0,0),(0,1,0)]
            self.NextNearestNeighbor[D][C]=[(0,0,0),(0,0,1)]

            # self.NextNearestNeighbor[A][B]=[(0,Ly-1,0),(0,0,Lz-1),(Lx-1,1,0),(Lx-1,0,1)]
            # self.NextNearestNeighbor[A][C]=[(Lx-1,0,0),(0,0,Lz-1),(1,Ly-1,0),(0,Ly-1,1)]
            # self.NextNearestNeighbor[A][D]=[(Lx-1,0,0),(0,Ly-1,0),(1,0,Lz-1),(0,1,Lz-1)]
            # self.NextNearestNeighbor[B][A]=[(1,Ly-1,0),(1,0,Lz-1),(0,0,1),(0,1,0)]
            # self.NextNearestNeighbor[B][C]=[(1,0,0),(1,0,Lz-1),(0,Ly-1,0),(0,Ly-1,1)]
            # self.NextNearestNeighbor[B][D]=[(1,Ly-1,0),(1,0,0),(0,0,Lz-1),(0,1,Lz-1)]
            # self.NextNearestNeighbor[C][A]=[(Lx-1,1,0),(0,1,Lz-1),(1,0,0),(0,0,1)]
            # self.NextNearestNeighbor[C][B]=[(0,1,0),(0,1,Lz-1),(Lx-1,0,0),(Lx-1,0,1)]
            # self.NextNearestNeighbor[C][D]=[(0,1,0),(Lx-1,1,0),(0,0,Lz-1),(1,0,Lz-1)]
            # self.NextNearestNeighbor[D][A]=[(Lx-1,0,1),(0,Ly-1,1),(1,0,0),(0,1,0)]
            # self.NextNearestNeighbor[D][B]=[(0,0,1),(0,Ly-1,1),(Lx-1,0,0),(Lx-1,1,0)]
            # self.NextNearestNeighbor[D][C]=[(0,0,1),(Lx-1,0,1),(0,Ly-1,0),(1,Ly-1,0),]

            # ########J3, only connected####################
            # self.NextNextNearestNeighbor[A][A]=[(1,0,0),(0,1,0),(0,0,1),(Lx-1,0,0),(0,Ly-1,0),(0,0,Lz-1)]
            # self.NextNextNearestNeighbor[B][B]=[(Lx-1,0,0),(Lx-1,1,0),(Lx-1,0,1),(1,0,0),(1,Ly-1,0),(1,0,Lz-1)]
            # self.NextNextNearestNeighbor[C][C]=[(0,Ly-1,0),(1,Ly-1,0),(0,Ly-1,1),(0,1,0),(Lx-1,1,0),(0,1,Lz-1)]
            # self.NextNextNearestNeighbor[D][D]=[(0,0,Lz-1),(1,0,Lz-1),(0,1,Lz-1),(0,0,1),(Lx-1,0,1),(0,Ly-1,1)]

            for i in range(6):
                for j in range(6):
                    for e in self.NearestNeighbor[i][j]:
                        #self.BareW.Data[:,i,:,j,self.__Map.CoordiIndex(e)]+= J1*SzSz;
                        self.BareW.Data[:,i,:,j,self.__Map.CoordiIndex(e)]+= J1*SS;
                    for e in self.NextNearestNeighbor[i][j]:
                        #self.BareW.Data[:,i,:,j,self.__Map.CoordiIndex(e)]+= J1*SzSz;
                        self.BareW.Data[:,i,:,j,self.__Map.CoordiIndex(e)]+= J2*SS;
                    # for e in self.NextNextNearestNeighbor[i][j]:
                    #     #self.BareW.Data[:,i,:,j,self.__Map.CoordiIndex(e)]+= J1*SzSz;
                    #     self.BareW.Data[:,i,:,j,self.__Map.CoordiIndex(e)]+= J3*SS;
            #S111S111=np.outer(Sx+Sy+Sz,Sx+Sy+Sz)
            #for i in range(4):
            #self.BareW.Data[:,i,:,i,0]+ = -10.0*J1*S111S111;

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

            for i in range(2):
                for j in range(2):
                    #J1 interaction A-->B, B-->A
                    for e in self.NearestNeighbor[i][j]:
                        self.BareW.Data[:,i,:,j,self.__Map.CoordiIndex(e)]+= J1*SS;
                    #J2 interaction A-->A, B-->B
                    for e in self.NextNearestNeighbor[i][j]:
                        self.BareW.Data[:,i,:,j,self.__Map.CoordiIndex(e)]+= J2*SS;
        else:
            Assert(False, "Lattice {0} has not been implemented yet!".format(LatName))


    def ToDict(self):
        points, lines=self.Lat.GetSitesList(SubLatIn=0)
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

