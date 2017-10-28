#!usr/bin/env python
import sys
import numpy as np
import parameter as para
from weight import UP,DOWN,IN,OUT
import weight, plot
from logger import *
import r_index 
# from scipy import weave
print "calculator"

def SimpleGG(G, _map):
    #half integer tin and tout
    GammaG=np.zeros([2, _map.Vol, _map.MaxTauBin, _map.MaxTauBin])+0.0*1j
    G.FFT("R", "T")
    Uspin=UP
    Gspin=UP
    sub=0
    r=0
    for t1 in range(_map.MaxTauBin):
        tg1=t1
        sign1=1
        G1=G.Data[Uspin,sub,Gspin,sub,r,tg1]*sign1
        for t2 in range(_map.MaxTauBin):
            tg2=-t2-1
            sign2=1
            if tg2<0:
                tg2+=_map.MaxTauBin
                sign2=-sign2
            G2=sign2*G.Data[Gspin,sub,Uspin,sub,r,tg2]
            GammaG[Gspin, r, t1, t2]=G1*G2 
    return GammaG

def GGW(GammaG,W,_map):
    W.FFT("R","T")
    #GammaG=SimpleGG(G,_map)
    OrderSign = -1
    spinUP=_map.Spin2Index(UP,UP)
    spinDOWN=_map.Spin2Index(DOWN,DOWN)
    spinUPDOWN=_map.Spin2Index(UP,DOWN)
    spinDOWNUP=_map.Spin2Index(DOWN,UP)
    sub=0
    r=0
    GGW=np.zeros([2, _map.Vol, _map.MaxTauBin, _map.MaxTauBin])+0.0*1j
    for t1 in range(_map.MaxTauBin):
        for t2 in range(_map.MaxTauBin):
            if t1==t2:
                print W.Data[spinUP,sub,spinUP,sub,0,abs(t1-t2)]
            GGW[UP,r, t1,t2]=OrderSign*GammaG[UP,r,t1,t2]*W.Data[spinUP,sub,spinUP,sub,0,abs(t1-t2)]
            GGW[DOWN,r, t1,t2]=OrderSign*GammaG[UP,r,t1,t2]*W.Data[spinDOWNUP,sub,spinUPDOWN,sub,0,abs(t1-t2)]
    return GGW

def AddTwoGToGammaG(GammaG, G, _map):
    #integer tin and tout
    G.FFT("R","T")
    spinUP=_map.Spin2Index(UP,UP)
    spinDOWN=_map.Spin2Index(DOWN,DOWN)
    sub=0
    r=0
    GGGammaG=np.zeros([2, _map.Vol, _map.MaxTauBin, _map.MaxTauBin])+0.0*1j
    for t in range(_map.MaxTauBin):
        for tout in range(_map.MaxTauBin):
            tgout=t-tout-1
            signout=1
            if tgout<0:
                tgout+=_map.MaxTauBin
                signout*=-1
            Gout=signout*G.Data[:,sub,:,sub,r,tgout]
            for tin in range(_map.MaxTauBin):
                GGGammaG[UP, r, t, tin]+=Gout[UP,UP]*GammaG[UP,r,tout,tin]
                GGGammaG[DOWN, r, t, tin]+=Gout[DOWN,DOWN]*GammaG[DOWN,r,tout,tin]
    GGGammaGNew=np.zeros([2, _map.Vol, _map.MaxTauBin, _map.MaxTauBin])+0.0*1j
    for t in range(_map.MaxTauBin):
        for tin in range(_map.MaxTauBin):
            tgin=tin-t 
            signin=1
            if tgin<0:
                tgin+=_map.MaxTauBin
                signin*=-1
            Gin=signin*G.Data[:,sub,:,sub,r,tgin]
            for tout in range(_map.MaxTauBin):
                GGGammaGNew[UP, r, tout, t]+=Gin[UP,UP]*GGGammaG[UP,r,tout,tin]
                GGGammaGNew[DOWN, r, tout, t]+=Gin[DOWN,DOWN]*GGGammaG[DOWN,r,tout,tin]
    return GGGammaGNew*_map.Beta**2/_map.MaxTauBin**2

def AddG_To_GammaG(GammaG, G, _map):
    #integer tin and tout
    G.FFT("R","T")
    FermiLoopSign=-1
    spinUP=_map.Spin2Index(UP,UP)
    spinDOWN=_map.Spin2Index(DOWN,DOWN)
    sub=0
    GGammaG = np.zeros([4, _map.Vol, _map.Vol, _map.MaxTauBin, _map.MaxTauBin])+0.0*1j
    for r in range(_map.Vol):
        for tout in range(_map.MaxTauBin):
            for tin in range(_map.MaxTauBin):
                tgout = tin - tout -1
                sign=1
                if tgout<0:
                    tgout+=_map.MaxTauBin
                    sign*=-1
                Gout = 0.5*sign*G.Data[:,sub,:,sub,0,tgout]

                tgout = tin - tout
                sign=1
                if tgout<0:
                    tgout+=_map.MaxTauBin
                    sign*=-1
                Gout += 0.5*sign*G.Data[:,sub,:,sub,0,tgout]

                ## UPUPUP
                GGammaG[0, r, r, tout, tin] = Gout[UP,UP]*GammaG[UP,r,tout,tin]
                ## UPDOWNUP
                GGammaG[1, r, r, tout, tin] = Gout[DOWN,DOWN]*GammaG[UP,r,tout,tin]
                ## DOWNUPDOWN
                GGammaG[2, r, r, tout, tin] = Gout[UP,UP]*GammaG[DOWN,r,tout,tin]
                ## DOWNDOWNDOWN
                GGammaG[3, r, r, tout, tin] = Gout[DOWN,DOWN]*GammaG[DOWN,r,tout,tin]
    return FermiLoopSign*GGammaG

def WWGammaW(GGammaG, W0, W, _map):
    sub = 0
    W.FFT("R","T")
    W0.FFT("R")
    UPUP=_map.Spin2Index(UP,UP)
    DOWNDOWN=_map.Spin2Index(DOWN,DOWN)
    UPDOWN=_map.Spin2Index(UP,DOWN)
    DOWNUP=_map.Spin2Index(DOWN,UP)

    print "calculating WGammaW..."
    WGammaW=np.zeros([6, _map.Vol, _map.Vol, _map.MaxTauBin, _map.MaxTauBin])+0.0*1j
    for r in range(_map.Vol):
        for rout in range(_map.Vol):
            dr_out = int(r_index.CoordiIndex(r, rout, _map))
            for t in range(_map.MaxTauBin):
                for tout in range(_map.MaxTauBin):

                    dt_out = t - tout -1
                    if dt_out<0:
                        dt_out+=_map.MaxTauBin

                    Wout = W.Data[:,sub,:,sub,dr_out,dt_out]*(_map.Beta/_map.MaxTauBin)
                    if t == tout:
                        Wout += W0.Data[:,sub,:,sub,dr_out]

                    ##For Delta W test
                    #if t!= tout:
                        #continue
                    #else:
                        #Wout = W0.Data[:,sub,:,sub,dr_out]

                    for rin in range(_map.Vol):
                        for tin in range(_map.MaxTauBin):
                            # UPUP UPUP
                            WGammaW[0, r, rin, t, tin]  += Wout[UPUP, UPUP] * GGammaG[0, rout, rin, tout, tin]

                            # DOWNDOWN DOWNDOWN
                            WGammaW[1, r, rin, t, tin]  += Wout[DOWNDOWN, DOWNDOWN] * GGammaG[3, rout, rin, tout, tin]

                            #UPUP DOWNDOWN 
                            WGammaW[2, r, rin, t, tin]  += Wout[UPUP, DOWNDOWN] * GGammaG[3, rout, rin, tout, tin]

                            #DOWNDOWN UPUP
                            WGammaW[3, r, rin, t, tin]  += Wout[DOWNDOWN, UPUP] * GGammaG[0, rout, rin, tout, tin]

                            # UPDOWN DOWNUP
                            WGammaW[4, r, rin, t, tin]  += Wout[UPDOWN, DOWNUP] * GGammaG[1, rout, rin, tout, tin]

                            # DOWNUP UPDOWN
                            WGammaW[5, r, rin, t, tin]  += Wout[DOWNUP, UPDOWN] * GGammaG[2, rout, rin, tout, tin]

    print "calculating WWGammaW..."
    WWGammaW = np.zeros([6, _map.Vol, _map.Vol, _map.MaxTauBin, _map.MaxTauBin]) + 0.0*1j
    for r in range(_map.Vol):
        for rin in range(_map.Vol):
            dr_in = int(r_index.CoordiIndex(rin, r, _map))
            for t in range(_map.MaxTauBin):
                for tin in range(_map.MaxTauBin):
                    dt_in = tin - t
                    if dt_in < 0:
                        dt_in +=_map.MaxTauBin

                    Win = W.Data[:,sub,:,sub,dr_in,dt_in]*(_map.Beta/_map.MaxTauBin)
                    if t == tin:
                        Win += W0.Data[:,sub,:,sub, dr_in]

                    ###For DeltaW test
                    #if t!= tin:
                        #continue
                    #else:
                        #Win = W0.Data[:,sub,:,sub,dr_in]

                    rout = r
                    for tout in range(_map.MaxTauBin):

                        ## UPUP UPUP
                        WWGammaW[0, rout, r, tout, t] += Win[UPUP, UPUP] * WGammaW[0, rout, rin, tout, tin]
                        WWGammaW[0, rout, r, tout, t] += Win[UPUP, DOWNDOWN] * WGammaW[2, rout, rin, tout, tin]
                        WWGammaW[0, r, rout, t, tout] += Win[UPUP, UPUP] * WGammaW[0, rout, rin, tout, tin]
                        WWGammaW[0, r, rout, t, tout] += Win[UPUP, DOWNDOWN] * WGammaW[2, rout, rin, tout, tin]

                        ## DOWNDOWN DOWNDOWN
                        WWGammaW[1, rout, r, tout, t] += Win[DOWNDOWN, UPUP] * WGammaW[3, rout, rin, tout, tin]
                        WWGammaW[1, rout, r, tout, t] += Win[DOWNDOWN, DOWNDOWN] * WGammaW[1, rout, rin, tout, tin]
                        WWGammaW[1, r, rout, t, tout] += Win[DOWNDOWN, UPUP] * WGammaW[3, rout, rin, tout, tin]
                        WWGammaW[1, r, rout, t, tout] += Win[DOWNDOWN, DOWNDOWN] * WGammaW[1, rout, rin, tout, tin]

                        ## UPUP DOWNDOWN
                        WWGammaW[2, rout, r, tout, t] += Win[DOWNDOWN, UPUP] * WGammaW[0, rout, rin, tout, tin]
                        WWGammaW[2, rout, r, tout, t] += Win[DOWNDOWN, DOWNDOWN] * WGammaW[2, rout, rin, tout, tin]
                        WWGammaW[2, r, rout, t, tout] += Win[UPUP, UPUP] * WGammaW[3, rout, rin, tout, tin]
                        WWGammaW[2, r, rout, t, tout] += Win[UPUP, DOWNDOWN] * WGammaW[1, rout, rin, tout, tin]

                        ## DOWNDOWN UPUP
                        WWGammaW[3, rout, r, tout, t] += Win[UPUP, UPUP] * WGammaW[3, rout, rin, tout, tin]
                        WWGammaW[3, rout, r, tout, t] += Win[UPUP, DOWNDOWN] * WGammaW[1, rout, rin, tout, tin]
                        WWGammaW[3, r, rout, t, tout] += Win[DOWNDOWN, UPUP] * WGammaW[0, rout, rin, tout, tin]
                        WWGammaW[3, r, rout, t, tout] += Win[DOWNDOWN, DOWNDOWN] * WGammaW[2, rout, rin, tout, tin]

                        ## UPDOWN DOWNUP
                        WWGammaW[4, rout, r, tout, t] += Win[DOWNUP, UPDOWN] * WGammaW[4, rout, rin, tout, tin]
                        WWGammaW[4, rout, r, tout, t] += Win[UPDOWN, DOWNUP] * WGammaW[5, rout, rin, tout, tin]

                        ## DOWNUP UPDOWN
                        WWGammaW[5, rout, r, tout, t] += Win[UPDOWN, DOWNUP] * WGammaW[5, rout, rin, tout, tin]
                        WWGammaW[5, r, rout, t, tout] += Win[DOWNUP, UPDOWN] * WGammaW[4, rout, rin, tout, tin]
    return -1.0*WWGammaW

def GammaWToGammaG(GammaW, G, _map):
    #integer tin and tout
    G.FFT("R","T")
    spinUP=_map.Spin2Index(UP,UP)
    spinDOWN=_map.Spin2Index(DOWN,DOWN)
    sub=0
    GGammaW = np.zeros([2, _map.Vol, _map.MaxTauBin, _map.MaxTauBin])+0.0*1j
    for r in range(_map.Vol):
        for tout in range(_map.MaxTauBin):
            for tin in range(_map.MaxTauBin):
                tgout = tout - tin -1
                sign=1
                if tgout<0:
                    tgout+=_map.MaxTauBin
                    sign*=-1
                Gout = 0.5*sign*G.Data[:,sub,:,sub,0,tgout]

                tgout = tout - tin
                sign=1
                if tgout<0:
                    tgout+=_map.MaxTauBin
                    sign*=-1
                Gout += 0.5*sign*G.Data[:,sub,:,sub,0,tgout]

                GGammaW[UP, r, tout, tin]+=Gout[UP,UP]*GammaW[0,r,r,tout,tin]
                GGammaW[UP, r, tout, tin]+=Gout[DOWN,DOWN]*GammaW[5,r,r,tout,tin]
                GGammaW[DOWN, r, tout, tin]+=Gout[DOWN,DOWN]*GammaW[1,r,r,tout,tin]
                GGammaW[DOWN, r, tout, tin]+=Gout[UP,UP]*GammaW[4,r,r,tout,tin]
    GammaG=AddTwoGToGammaG(GGammaW, G, _map)
    return GammaG

def shift(r, L):
    if r<0:
        r+=L
    elif r>=L:
        r-=L
    return r

def GammaG_FirstOrder(GammaG, G, W0, _map):
    #OrderSign=-1, FermiLoopSign=-1, therefore TotalSign=1
    #integer tin and tout
    GG=SimpleGG(G,_map)
    GammaGNew=np.zeros([2, _map.Vol, _map.MaxTauBin, _map.MaxTauBin])+0.0*1j
    spinUP=_map.Spin2Index(UP,UP)
    spinDOWN=_map.Spin2Index(DOWN,DOWN)

    # # # print "GammaG[UP,UP]=\n", GammaG[UP,0,:,-1]
    W0.FFT("R")
    sub=0
    r=0
    Neighbors=[]
    Lx, Ly=_map.L
    for Gx in range(Lx):
        for dx in [-1,0,1]:
            x=shift(Gx+dx, Lx)
            #if x<0:
                #x+=Lx
            #elif x>=Lx:
                #x-=Lx
            for Gy in range(Ly):
                for dy in [-1,0,1]:
                    y=shift(Gy+dy, Ly)
                    #if y<0:
                        #y+=Ly
                    #elif y>=Ly:
                        #y-=Ly

                    dx_shift=shift(dx, Lx)
                    dy_shift=shift(dy, Ly)
                    #if dx<0
                        #dx_shift+=Lx
                    #if dy<0:
                        #dy_shift+=Ly

                    i=_map.CoordiIndex([Gx,Gy])
                    j=_map.CoordiIndex([x,y])
                    k=_map.CoordiIndex([dx_shift,dy_shift])
                    Neighbors.append([i,j,W0.Data[:,sub,:,sub,k]])
                    # print W0.Data[spinUP,sub,spinUP,sub,k], dx, dy
                    # print W0.Data[:,sub,:,sub,k], dx, dy

    for t3 in range(_map.MaxTauBin):
        for tin in range(_map.MaxTauBin):
            dtin=t3-tin
            sign=1
            if dtin<0:
                dtin+=_map.MaxTauBin
                sign*=-1
            G1=sign*G.Data[UP,sub,UP,sub,0,dtin]

            for tout in range(_map.MaxTauBin):
                dtout=tout-t3-1
                sign=1
                if dtout<0:
                    dtout+=_map.MaxTauBin
                    sign*=-1
                G2=sign*G.Data[UP,sub,UP,sub,0,dtout]
                GG=G1*G2
                
                for r1,r2,V in Neighbors:
                    GammaGNew[UP,r1,tout,tin]+=GG*(V[spinUP,spinUP]*GammaG[UP,r2,t3,t3]+V[spinUP,spinDOWN]*GammaG[DOWN,r2,t3,t3])
                    GammaGNew[DOWN,r1,tout,tin]+=GG*(V[spinDOWN,spinUP]*GammaG[UP,r2,t3,t3]+V[spinDOWN,spinDOWN]*GammaG[DOWN,r2,t3,t3])

    GammaGNew*=_map.Beta/_map.MaxTauBin
    return GammaGNew

def SigmaSmoothT_FirstOrder(G, W, map):
    '''Fock diagram, assume Spin Conservation'''
    OrderSign=-1
    Sigma=weight.Weight("SmoothT", map, "TwoSpins", "AntiSymmetric", "R", "T")
    TauRange = range(map.MaxTauBin)
    G.FFT("R", "T")
    W.FFT("R", "T")
    for spin1 in range(2):
        for spin2 in range(2):
            spinW = (map.Spin2Index(spin1,spin2),map.Spin2Index(spin2,spin1))
            spinG = (spin2, spin2)
            spinSigma = (spin1, spin1)
            Sigma.Data[spinSigma[IN], :, spinSigma[OUT], :, :]  \
                    += OrderSign*G.Data[spinG[IN], :, spinG[OUT], :, :]\
                    *W.Data[spinW[IN], :, spinW[OUT], :, :]
    return Sigma

def SigmaDeltaT_FirstOrder(G, W0, map):
    '''Hatree-Fock diagram, assume Spin Conservation'''
    ########Fock Diagram
    OrderSign=-1
    AntiSymmetricFactor=-1
    SigmaDeltaT=weight.Weight("DeltaT", map, "TwoSpins", "AntiSymmetric", "R")
    G.FFT("R", "T")
    W0.FFT("R")
    for spin1 in range(2):
        for spin2 in range(2):
            #############G(tau==-0)
            spinW = (map.Spin2Index(spin1,spin2),map.Spin2Index(spin2,spin1))
            spinG = (spin2, spin2)
            spinSigma = (spin1, spin1)
            SigmaDeltaT.Data[spinSigma[IN], :, spinSigma[OUT], :, :]+= OrderSign \
                    *AntiSymmetricFactor*(1.5*G.Data[spinG[IN], :, spinG[OUT], :, :, -1]\
                    -0.5*G.Data[spinG[IN], :,spinG[OUT], :, :,-2])\
                    *W0.Data[spinW[IN], :, spinW[OUT], :, :]

    ########Hatree Diagram, or bubble diagram
    FermiLoopSign=-1
    for sp1 in range(2):
        for sp2 in range(2):
            spinW = (map.Spin2Index(sp1,sp1), map.Spin2Index(sp2,sp2))
            for sub1 in range(map.NSublat):
                for sub2 in range(map.NSublat):
                    for r in range(map.Vol):
                        G1 = 1.5*G.Data[sp2,sub2,sp2,sub2,0,-1]-0.5*G.Data[sp2,sub2,sp2,sub2,0,-2]
                        SigmaDeltaT.Data[sp1, sub1, sp1, sub1, 0]+= OrderSign*FermiLoopSign \
                            *AntiSymmetricFactor*G1*W0.Data[spinW[IN], sub1, spinW[OUT], sub2, r]
    return SigmaDeltaT


def Polar_FirstOrder(G, map):
    OrderSign=-1
    FermiLoopSign=-1
    AntiSymmetricFactor=-1
    Polar=weight.Weight("SmoothT", map, "FourSpins","Symmetric", "R","T")
    G.FFT("R","T")
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
                    map.Spin2Index(*spinPolar[OUT]),subB,:,:]+=OrderSign*FermiLoopSign \
                    *AntiSymmetricFactor*G.Data[spinG1[IN], subB, spinG1[OUT], subA, :, ::-1]  \
                    *G.Data[spinG2[IN], subA, spinG2[OUT], subB, :, :]
    return Polar

def W_FirstOrder(W0, Polar, map):
    W=weight.Weight("SmoothT", map, "FourSpins", "Symmetric", "K","W")
    W0.FFT("K")
    Polar.FFT("K","W")
    Beta=map.Beta
    TauRange = range(map.MaxTauBin)
    SubRange=range(W.NSublat)
    SubList=[(a,b,c,d) for a in SubRange for b in SubRange for c in SubRange for d in SubRange]
    SpinList=[(Wtuple,Polartuple) for Wtuple in map.GetConservedSpinTuple("FourSpins") \
                                  for Polartuple in map.GetConservedSpinTuple("FourSpins") \
                                  if map.IsConserved(4, (Wtuple[IN], Polartuple[IN]))] 
    #make sure spin conservation on W0

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
    return W

def Calculate_Denorminator(W0, Polar, map):
    """require W0 and Polar to be in k, omega space"""
    W0.FFT("K")
    Polar.FFT("K","W")
    Beta=map.Beta
    NSpin, NSub=W0.NSpin, W0.NSublat
    JP=np.einsum("ijklv,klmnvt->ijmnvt",W0.Data, Polar.Data)
    #JP shape: NSpin,NSub,NSpin,NSub,Vol,Tau
    Temp=JP*Beta/map.MaxTauBin
    for tau in range(map.MaxTauBin):
        Temp[...,tau] *= np.cos(tau*np.pi/map.MaxTauBin)
    I=np.eye(NSpin*NSub).reshape([NSpin,NSub,NSpin,NSub])
    return I[...,np.newaxis,np.newaxis]-Temp, JP

def W_Dyson(W0, Polar, map, Lat):
    W=weight.Weight("SmoothT", map, "FourSpins", "Symmetric", "K","W")
    ChiTensor=weight.Weight("SmoothT", map, "FourSpins", "Symmetric", "K","W")

    Polar.FFT("K","W")
    Denorm,JP=Calculate_Denorminator(W0, Polar, map)

    W0.FFT("K")
    JPJ=np.einsum("ijklvt,klmnv->ijmnvt", JP, W0.Data)
    lu_piv,Determ=weight.LUFactor(Denorm)
    Check_Denorminator(Denorm, Determ,map)
    W.LUSolve(lu_piv, JPJ)
    ChiTensor.LUSolve(lu_piv, -Polar.Data)
    return W, ChiTensor, Determ

def G_Dyson(G0, SigmaDeltaT, Sigma, map):
    Beta=map.Beta
    G=weight.Weight("SmoothT", map, "TwoSpins", "AntiSymmetric", "K","W")
    G0.FFT("K", "W")
    SigmaDeltaT.FFT("K")
    Sigma.FFT("K", "W")

    NSpin, NSub=G.NSpin, G.NSublat

    G0SigmaDeltaT=np.einsum("ijklvt,klmnv->ijmnvt",G0.Data, SigmaDeltaT.Data)
    G0Sigma=np.einsum("ijklvt,klmnvt->ijmnvt",G0.Data, Sigma.Data)

    ####correction term
    for tau in range(map.MaxTauBin):
        G0SigmaDeltaT[...,tau]*= np.cos(np.pi*map.IndexToTau(tau)/Beta)

    GS  = Beta/map.MaxTauBin*(Beta/map.MaxTauBin*G0Sigma+G0SigmaDeltaT) 
    #GS shape: NSpin,NSub,NSpin,NSub,Vol,Tau

    I=np.eye(NSpin*NSub).reshape([NSpin,NSub,NSpin,NSub])
    Denorm=I[...,np.newaxis,np.newaxis]-GS
    lu_piv,Determ=weight.LUFactor(Denorm)
    Check_Denorminator(Denorm, Determ,map)
    G.LUSolve(lu_piv, G0.Data);
    return G

def Add_ChiTensor_ZerothOrder(ChiTensor, G, map):
    """add G(tau=0^-, k=0)G(tau=0^-, k=0) to ChiTensor
       This diagram has +1 sign since two fermi loop contribute (-1)^2
    """
    G.FFT("R", "T")
    ChiTensor.FFT("R", "T")
    NSub=map.NSublat
    SpList=[(a,b,c,d) for a in range(2) for b in range(2) for c in range(2) for d in range(2)]
    SubList=[(a,b) for a in range(NSub) for b in range(NSub)]
    ChiTensorInUnitCell=np.empty((ChiTensor.NSpin,NSub,ChiTensor.NSpin,NSub), dtype=complex)
    for spIn1, spIn2, spOut1, spOut2 in SpList:
        for subIn, subOut in SubList:
            for tau in range(map.MaxTauBin):
                ChiSpIn=map.Spin2Index(spIn1,spIn2)
                ChiSpOut=map.Spin2Index(spOut1,spOut2)
                G1 = 1.5*G.Data[spIn1,subIn,spIn2,subIn,0,-1]-0.5*G.Data[spIn1,subIn,spIn2,subIn,0,-2]
                G2 = 1.5*G.Data[spOut1,subOut,spOut2,subOut,0,-1]-0.5*G.Data[spOut1,subOut,spOut2,subOut,0,-2]
                ChiTensorInUnitCell[ChiSpIn,subIn,ChiSpOut,subOut]= G1 * G2
                #if spIn1==UP and spIn2==UP and spOut1==DOWN and spOut2==DOWN:
                    #if abs(G1)>0.1 and abs(G2)>0.1:
                        #print G1,G2, G1*G2
    ChiTensor.Data+=ChiTensorInUnitCell[...,np.newaxis,np.newaxis]
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
    Chi=weight.Weight("SmoothT", map, "NoSpin", "Symmetric", ChiTensor.SpaceDomain, ChiTensor.TimeDomain)
    Chi.Data=0.0
    #SS=[SxSx/4.0, SySy/4.0, SzSz/4.0]
    #for i in range(3):
        #temp=np.einsum("ik, kminvt->mnvt", SS[i], ChiTensor.Data)
        #Chi.Data+=temp.reshape([1, NSublat, 1, NSublat, map.Vol, map.MaxTauBin]) 

    SS=[SzSz/4.0]
    for i in range(len(SS)):
        temp=np.einsum("ik, kminvt->mnvt", SS[i], ChiTensor.Data)
        Chi.Data+=temp.reshape([1, NSublat, 1, NSublat, map.Vol, map.MaxTauBin]) 
    return Chi

class DenorminatorTouchZero(Exception):
    def __init__(self, value, pos, freq):
       self.value = value
       self.position=pos
       self.frequency=freq
    def __str__(self):
        return "Denorminator touch zero, minmum {0} is at K={1} and Omega={2}".format(self.value, self.position, self.frequency)

def Check_Denorminator(Denorm, Determ, map):
    pos=np.where(Determ==Determ.min())
    x,t=pos[0][0], pos[1][0]
    log.info("The minmum {0} is at K={1} and Omega={2}".format(Determ.min(), map.IndexToCoordi(x), t))
    SpSub,Vol,Time=Denorm.shape[0]*Denorm.shape[1], Denorm.shape[-2], Denorm.shape[-1]
    Denorm=Denorm.reshape([SpSub,SpSub,Vol,Time])
    log.info("The 1/linalg.cond is {0}".format(1.0/np.linalg.cond(Denorm[...,x,t])))
    if Determ.min().real<0.0 and Determ.min().imag<1.0e-4:
        raise DenorminatorTouchZero(Determ.min(), map.IndexToCoordi(x), t)
