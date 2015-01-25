#!usr/bin/env python
import sys
import numpy as np
import parameter as para
from weight import UP,DOWN,IN,OUT
import weight
from logger import *

def SigmaSmoothT_FirstOrder(G, W, map):
    '''Fock diagram'''
    Sigma=weight.Weight("SmoothT", map, "TwoSpins", "AntiSymmetric")

    TauRange = range(map.MaxTauBin)

    for spin1 in range(2):
        for spin2 in range(2):
            spinW = (map.Spin2Index(spin1,spin2),map.Spin2Index(spin2,spin1))
            spinG = (spin2, spin2)
            spinSigma = (spin1, spin1)
            Sigma.Data[spinSigma[IN], :, spinSigma[OUT], :, :]  \
                    -= G.Data[spinG[IN], :, spinG[OUT], :, :]\
                    *W.Data[spinW[IN], :, spinW[OUT], :, :]
    return Sigma

def SigmaDeltaT_FirstOrder(G, W0, map):
    '''Hatree-Fock diagram'''
    ########Fock Diagram
    Sigma0=weight.Weight("DeltaT", map, "TwoSpins", "AntiSymmetric")
    for spin1 in range(2):
        for spin2 in range(2):
            #############G(tau==-0)
            spinW = (map.Spin2Index(spin1,spin2),map.Spin2Index(spin2,spin1))
            spinG = (spin2, spin2)
            spinSigma = (spin1, spin1)
            Sigma0.Data[spinSigma[IN], :, spinSigma[OUT], :, :]  \
                    -= G.Data[spinG[IN], :, spinG[OUT], :, :, -1]\
                    *W0.Data[spinW[IN], :, spinW[OUT], :, :]
    ########Hatree Diagram
    for spin1 in range(2):
        for spin2 in range(2):
            spinW = (map.Spin2Index(spin1,spin1), map.Spin2Index(spin2,spin2))
            spinG = (spin2,spin2)
            spinSigma = (spin1, spin1)
            for r in range(map.Vol):
                Sigma0.Data[spinSigma[IN], :, spinSigma[OUT], :, 0] \
                        -= -G.Data[spinG[IN], :, spinG[OUT], :, 0, -1] \
                        *W0.Data[spinW[IN], :, spinW[OUT], :, r]
    return Sigma0


def Polar_FirstOrder(G, map):
    Polar=weight.Weight("SmoothT", map, "FourSpins","Symmetric")
    NSublat = map.NSublat
    NSpin = G.NSpin
    SubList=[(a,b) for a in range(NSublat) for b in range(NSublat)]
    SpList=[(a,b) for a in range(NSpin) for b in range(NSpin)]
    for spin1,spin2 in SpList:
        spinPolar = ((spin1,spin2),(spin2,spin1))
        spinG1 = (spin1, spin1)
        spinG2 = (spin2, spin2)
        for subA,subB in SubList:
            Polar.Data[map.Spin2Index(*spinPolar[IN]),subA, \
                    map.Spin2Index(*spinPolar[OUT]),subB,:,:]\
                    = (-1.0)*G.Data[spinG1[IN], subB, spinG1[OUT], subA, :, ::-1]  \
                    *G.Data[spinG2[IN], subA, spinG2[OUT], subB, :, :]
    return Polar

def W_FirstOrder(W0, Polar, map):
    W=weight.Weight("SmoothT", map, "FourSpins", "Symmetric")
    Beta=map.Beta
    TauRange = range(map.MaxTauBin)
    SubRange=range(W.NSublat)
    SubList=[(a,b,c,d) for a in SubRange for b in SubRange for c in SubRange for d in SubRange]
    SpinList=[(Wtuple,Polartuple) for Wtuple in map.GetConservedSpinTuple("FourSpins") \
                                  for Polartuple in map.GetConservedSpinTuple("FourSpins") \
                                  if map.IsConserved(4, (Wtuple[IN], Polartuple[IN]))] 
    #make sure spin conservation on W0

    W0.FFT(1, "Space")
    Polar.FFT(1, "Space")

    for spWt,spPolart in SpinList:
        spW0L=(map.Spin2Index(*spWt[IN]), map.Spin2Index(*spPolart[IN]))
        spW0R=(map.Spin2Index(*spPolart[OUT]), map.Spin2Index(*spWt[IN]))
        spW = (map.Spin2Index(*spWt[IN]), map.Spin2Index(*spWt[OUT]))
        spPolar = (map.Spin2Index(*spPolart[IN]), map.Spin2Index(*spPolart[OUT]))
        for e in SubList:
            for tau in TauRange:
                W.Data[spW[IN],e[0],spW[OUT],e[3],:,tau]+=\
                        W0.Data[spW0L[IN],e[0],spW0L[OUT],e[1],:] \
                        *Polar.Data[spPolar[IN],e[1],spPolar[OUT],e[2],:,tau]\
                        *W0.Data[spW0R[IN],e[2],spW0R[OUT],e[3],:]
    
    W0.FFT(-1, "Space")
    Polar.FFT(-1, "Space")
    W.FFT(-1, "Space")
    return W

def Calculate_Denorminator(W0, Polar, map):
    """require W0 and Polar to be in k, omega space"""
    Beta=map.Beta
    NSpin, NSub=W0.NSpin, W0.NSublat
    JP=np.einsum("ijklv,klmnvt->ijmnvt",W0.Data, Polar.Data)
    #JP shape: NSpin,NSub,NSpin,NSub,Vol,Tau
    JP *= Beta/map.MaxTauBin
    for tau in range(map.MaxTauBin):
        JP[...,tau] *= np.cos(tau*np.pi/map.MaxTauBin)
    I=np.eye(NSpin*NSub).reshape([NSpin,NSub,NSpin,NSub])
    return I[...,np.newaxis,np.newaxis]-JP

def W_Dyson(W0, Polar, map):
    W=weight.Weight("SmoothT", map, "FourSpins", "Symmetric")
    W0.FFT(1, "Space")
    Polar.FFT(1, "Space", "Time")

    JP=np.einsum("ijklv,klmnvt->ijmnvt",W0.Data, Polar.Data)
    JPJ=np.einsum("ijklvt,klmnv->ijmnvt", JP, W0.Data)
    W.Data=Calculate_Denorminator(W0, Polar, map)
    W.Inverse();
    W.Data = np.einsum('ijklvt,klmnvt->ijmnvt', W.Data,JPJ)

    W.FFT(-1, "Space", "Time")
    Polar.FFT(-1, "Space", "Time")
    W0.FFT(-1, "Space")
    return W

def G_Dyson(G0, Sigma0, Sigma, map):
    Beta=map.Beta
    G=weight.Weight("SmoothT", map, "TwoSpins", "AntiSymmetric")
    G0.FFT(1, "Space", "Time")
    Sigma0.FFT(1, "Space")
    Sigma.FFT(1, "Space", "Time")

    NSpin, NSub=G.NSpin, G.NSublat

    G0Sigma0=np.einsum("ijklvt,klmnv->ijmnvt",G0.Data, Sigma0.Data)
    G0Sigma=np.einsum("ijklvt,klmnvt->ijmnvt",G0.Data, Sigma.Data)

    ####correction term
    for tau in range(map.MaxTauBin):
        G0Sigma0[...,tau]*= np.cos(np.pi*map.IndexToTau(tau)/Beta)

    #GS  = Beta/map.MaxTauBin*(Beta/map.MaxTauBin*G0Sigma + G0Sigma0)
    GS  = Beta/map.MaxTauBin*(Beta/map.MaxTauBin*G0Sigma) 
    #GS shape: NSpin,NSub,NSpin,NSub,Vol,Tau

    I=np.eye(NSpin*NSub).reshape([NSpin,NSub,NSpin,NSub])
    G.Data=I[...,np.newaxis,np.newaxis]-GS
    G.Inverse();
    G.Data=np.einsum('ijklvt,klmnvt->ijmnvt', G.Data,G0.Data)

    G.FFT(-1, "Space", "Time")
    Sigma0.FFT(-1, "Space")
    Sigma.FFT(-1, "Space", "Time")
    G0.FFT(-1, "Space","Time")
    return G


def Calculate_ChiTensor(W0, Polar, map):
    ChiTensor=weight.Weight("SmoothT", map, "FourSpins", "Symmetric")

    W0.FFT(1, "Space")
    Polar.FFT(1, "Space", "Time")

    ChiTensor.Data=Calculate_Denorminator(W0, Polar, map)
    ChiTensor.Inverse();
    ChiTensor.Data = -np.einsum('ijklvt,klmnvt->ijmnvt', Polar.Data, ChiTensor.Data)

    ChiTensor.FFT(-1, "Space", "Time")
    Polar.FFT(-1, "Space", "Time")
    W0.FFT(-1, "Space")
    return ChiTensor

def Calculate_Chi(ChiTensor, map):
    NSpin, NSublat=ChiTensor.NSpin, ChiTensor.NSublat
    SxSx=np.zeros((NSpin,NSpin))
    SySy=np.zeros((NSpin,NSpin))
    SzSz=np.zeros((NSpin,NSpin))
    UU=map.Spin2Index(UP,UP) 
    UD=map.Spin2Index(UP,DOWN) 
    DU=map.Spin2Index(DOWN, UP) 
    DD=map.Spin2Index(DOWN, DOWN) 
    SxSx[UD, UD]=SxSx[DU, DU]=1
    SxSx[UD, DU]=SxSx[DU, UD]=1
    SySy[UD, UD]= SySy[DU, DU]=-1
    SySy[UD, DU]= SySy[DU, UD]=1
    SzSz[UU, UU]= SzSz[DD, DD]=1
    SzSz[UU, DD]= SzSz[DD, UU]=-1
    Chi=weight.Weight("SmoothT", map, "NoSpin", "Symmetric")
    Chi_ss=[Chi.Copy(), Chi.Copy(), Chi.Copy()]
    SS=[SxSx/4.0, SySy/4.0, SzSz/4.0]
    for i in range(3):
        temp=np.einsum("ik, kminvt->mnvt", SS[i], ChiTensor.Data)
        Chi_ss[i].Data=temp.reshape([1, NSublat, 1, NSublat, map.Vol, map.MaxTauBin]) 
    Chi.Data=Chi_ss[0].Data+Chi_ss[1].Data+Chi_ss[2].Data
    return Chi, Chi_ss

def Check_Denorminator(W0, Polar, map):
    """return tuple ((position, smallest 1-JP determinant), 1-JP in omega,k domain)"""
    log.info("Check Denorminator...")
    W0.FFT(1, "Space")
    Polar.FFT(1, "Space", "Time")
    NSpin, NSub=Polar.NSpin, Polar.NSublat
    Denorm=Calculate_Denorminator(W0, Polar, map)
    Denorm=Denorm.reshape([NSpin*NSub, NSpin*NSub]+Polar.Shape[Polar.VOLDIM:])
    Determ=np.zeros((map.Vol,map.MaxTauBin))*1j
    for x in range(map.Vol):
        for t in range(map.MaxTauBin):
            Determ[x,t]=np.linalg.det(Denorm[:,:,x,t])
    pos=np.where(Determ==Determ.min())
    x,t=pos[0][0], pos[1][0]
    log.info("The minmum {0} is at K={1} and Omega={2}".format(Determ.min(), map.IndexToCoordi(x), t))
    W0.FFT(-1, "Space")
    Polar.FFT(-1, "Space", "Time")
    return (x, t, Determ.min()), Determ
