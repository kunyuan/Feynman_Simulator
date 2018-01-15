#!usr/bin/env python
import sys
import numpy as np
import parameter as para
from weight import UP,DOWN,IN,OUT
import weight, plot
from logger import *
import r_index 
import calculator as calc
import gc

def SimpleGG(G, _map):
    #half integer tin and tout
    GGGammaG=np.zeros([2, _map.Vol, _map.MaxTauBin, _map.MaxTauBin])+0.0*1j
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
            GGGammaG[Gspin, r, t1, t2]=G1*G2 
    return GGGammaG

def AddW_To_GGGammaG(GGGammaG,W,_map):
    W.FFT("R","T")
    OrderSign = -1
    spinUP=_map.Spin2Index(UP,UP)
    spinDOWN=_map.Spin2Index(DOWN,DOWN)
    spinUPDOWN=_map.Spin2Index(UP,DOWN)
    spinDOWNUP=_map.Spin2Index(DOWN,UP)
    sub=0
    GammaG=np.zeros([2, _map.Vol, _map.MaxTauBin, _map.MaxTauBin])+0.0*1j
    Wshift=weight.Weight("SmoothT", _map, "FourSpins", "Symmetric", "R","T")
    for t in range(_map.MaxTauBin):
        t1=t-1
        if t1<0:
            t1+=_map.MaxTauBin
        Wshift.Data[:,:,:,:,:,t]=0.5*(W.Data[:,:,:,:,:,t1]+W.Data[:,:,:,:,:,t])

    for t1 in range(_map.MaxTauBin):
        for t2 in range(_map.MaxTauBin):
            dt = t1 - t2
            if dt <0:
                dt = dt + _map.MaxTauBin

            GammaG[UP, :, t1,t2] = GGGammaG[UP,:,t1,t2]*Wshift.Data[spinUP,sub,spinUP,sub,0,dt]
            GammaG[UP, :, t1,t2] += GGGammaG[DOWN,:,t1,t2]*Wshift.Data[spinUPDOWN,sub,spinDOWNUP,sub,0,dt]

            GammaG[DOWN, :, t1,t2] = GGGammaG[UP,:,t1,t2]*Wshift.Data[spinDOWNUP,sub,spinUPDOWN,sub,0,dt]
            GammaG[DOWN, :, t1,t2] += GGGammaG[DOWN,:,t1,t2]*Wshift.Data[spinDOWN,sub,spinDOWN,sub,0,dt]
    return GammaG*OrderSign

def AddTwoG_To_GammaG(GammaG, G, _map):
    #integer tin and tout
    G.FFT("R","T")
    spinUP=_map.Spin2Index(UP,UP)
    spinDOWN=_map.Spin2Index(DOWN,DOWN)
    sub=0
    r=0

    GGammaG=np.zeros([2, _map.Vol, _map.MaxTauBin, _map.MaxTauBin])+0.0*1j
    for t in range(_map.MaxTauBin):
        for tout in range(_map.MaxTauBin):
            tgout=t-tout-1
            signout=1
            if tgout<0:
                tgout+=_map.MaxTauBin
                signout*=-1
            Gout=0.5*signout*G.Data[:,sub,:,sub,r,tgout]
            tgout=t-tout
            signout=1
            if tgout<0:
                tgout+=_map.MaxTauBin
                signout*=-1
            Gout+=0.5*signout*G.Data[:,sub,:,sub,r,tgout]
            for tin in range(_map.MaxTauBin):
                GGammaG[UP, :, t, tin]+=Gout[UP,UP]*GammaG[UP,:,tout,tin]
                GGammaG[DOWN, :, t, tin]+=Gout[DOWN,DOWN]*GammaG[DOWN,:,tout,tin]

    GGGammaG=np.zeros([2, _map.Vol, _map.MaxTauBin, _map.MaxTauBin])+0.0*1j
    for t in range(_map.MaxTauBin):
        for tin in range(_map.MaxTauBin):
            tgin=tin-t 
            signin=1
            if tgin<0:
                tgin+=_map.MaxTauBin
                signin*=-1
            Gin=0.5*signin*G.Data[:,sub,:,sub,r,tgin]
            tgin=tin-t-1 
            signin=1
            if tgin<0:
                tgin+=_map.MaxTauBin
                signin*=-1
            Gin+=0.5*signin*G.Data[:,sub,:,sub,r,tgin]
            for tout in range(_map.MaxTauBin):
                GGGammaG[UP, :, tout, t]+=Gin[UP,UP]*GGammaG[UP,:,tout,tin]
                GGGammaG[DOWN, :, tout, t]+=Gin[DOWN,DOWN]*GGammaG[DOWN,:,tout,tin]
    return GGGammaG*_map.Beta**2/_map.MaxTauBin**2

def AddG_To_GGGammaG(GGGammaG, G, _map):
    #integer tin and tout
    G.FFT("R","T")
    FermiLoopSign=-1
    spinUP=_map.Spin2Index(UP,UP)
    spinDOWN=_map.Spin2Index(DOWN,DOWN)
    sub=0
    GammaW = np.zeros([2, _map.Vol, _map.Vol, _map.MaxTauBin, _map.MaxTauBin], dtype=np.complex)
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

            for r in range(_map.Vol):
                ## UP UP UP UP
                #GammaW[0, r, r, tout, tin] = Gout[UP,UP]*GGGammaG[UP,r,tout,tin]

                ## DOWN DOWN DOWN DOWN
                #GammaW[1, r, r, tout, tin] = Gout[DOWN,DOWN]*GGGammaG[DOWN,r,tout,tin]

                ## in:UP DOWN out:DOWN UP
                #GammaW[5, r, r, tout, tin] = Gout[UP, UP]*GGGammaG[DOWN,r,tout,tin]

                ## in:DOWN UP out:UP DOWN 
                #GammaW[4, r, r, tout, tin] = Gout[DOWN,DOWN]*GGGammaG[UP,r,tout,tin]

                # UP UP UP UP+ DOWN DOWN DOWN DOWN
                GammaW[0, r, r, tout, tin] = Gout[UP,UP]*GGGammaG[UP,r,tout,tin]+Gout[DOWN,DOWN]*GGGammaG[DOWN,r,tout,tin]

                # in:DOWN UP out:UP DOWN 
                GammaW[1, r, r, tout, tin] = Gout[DOWN,DOWN]*GGGammaG[UP,r,tout,tin]

            ### reverse
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

            for r in range(_map.Vol):
                ## UP UP UP UP
                #GammaW[0, r, r, tout, tin] += Gout[UP,UP]*GGGammaG[UP,r,tin,tout]

                ## DOWN DOWN DOWN DOWN
                #GammaW[1, r, r, tout, tin] += Gout[DOWN,DOWN]*GGGammaG[DOWN,r,tin,tout]

                ## in:UP DOWN out:DOWN UP
                #GammaW[5, r, r, tout, tin] += Gout[DOWN, DOWN]*GGGammaG[UP,r,tin,tout]

                ## in:DOWN UP out:UP DOWN 
                #GammaW[4, r, r, tout, tin] += Gout[UP,UP]*GGGammaG[DOWN,r,tin,tout]

                ## UP UP UP UP + DOWN DOWN DOWN DOWN
                GammaW[0, r, r, tout, tin] += Gout[UP,UP]*GGGammaG[UP,r,tin,tout]+Gout[DOWN,DOWN]*GGGammaG[DOWN,r,tin,tout]

                ## in:DOWN UP out:UP DOWN 
                GammaW[1, r, r, tout, tin] += Gout[UP,UP]*GGGammaG[DOWN,r,tin,tout]
    return FermiLoopSign*GammaW

def GenerateSpinIndex(_map):
    UPUP=_map.Spin2Index(UP,UP)
    DOWNDOWN=_map.Spin2Index(DOWN,DOWN)
    UPDOWN=_map.Spin2Index(UP,DOWN)
    DOWNUP=_map.Spin2Index(DOWN,UP)

    spinindex=np.array([UP,DOWN])
    spin2index=np.array([UPUP,DOWNDOWN, UPDOWN, DOWNUP])
    return (spinindex, spin2index)

def FFTGammaW(GammaW, _map, BackForth):
    OldShape=GammaW.shape
    NewShape=(_map.L[0], _map.L[1], _map.L[0], _map.L[1], _map.MaxTauBin, _map.MaxTauBin)
    GammaW=GammaW.reshape(NewShape)
    if BackForth==1:
        GammaW=np.fft.fftn(GammaW, axes=[0,1,2,3,4,5]) 
    elif BackForth==-1:
        GammaW=np.fft.ifftn(GammaW, axes=[0,1,2,3,4,5]) 
    GammaW=GammaW.reshape(OldShape)
    return GammaW

def FFTWshift(Wshift, _map, BackForth):
    OldShape=Wshift.shape
    NewShape=(4,4, _map.L[0], _map.L[1], _map.MaxTauBin)
    Wshift=Wshift.reshape(NewShape)
    if BackForth==1:
        Wshift=np.fft.fftn(Wshift, axes=[2,3,4]) 
    elif BackForth==-1:
        Wshift=np.fft.ifftn(Wshift, axes=[2,3,4]) 
    Wshift=Wshift.reshape(OldShape)
    return Wshift

def AddTwoW_To_GammaW(GammaW, W0, W, _map):
    # import gamma3
    sub = 0
    UPUP=_map.Spin2Index(UP,UP)
    DOWNDOWN=_map.Spin2Index(DOWN,DOWN)
    UPDOWN=_map.Spin2Index(UP,DOWN)
    DOWNUP=_map.Spin2Index(DOWN,UP)

    spinindex, spin2index=GenerateSpinIndex(_map)

    W.FFT("R","T")
    W0.FFT("R", "T")
    Wshift=weight.Weight("SmoothT", _map, "FourSpins", "Symmetric", "R","T")
    for t in range(_map.MaxTauBin):
        t1=t-1
        if t1<0:
            t1+=_map.MaxTauBin
        Wshift.Data[:,:,:,:,:,t]=0.5*(W.Data[:,:,:,:,:,t1]+W.Data[:,:,:,:,:,t])
    Wtot=np.array(Wshift.Data[:,0,:,0,:,:])*_map.Beta/_map.MaxTauBin
    Wtot[:,:,:,0]+=W0.Data[:,0,:,0,:]
    Wtot=FFTWshift(Wtot, _map, 1)
    #Wtot shape: spin1, spin2, dr, dt

    WWGammaW=np.zeros([2, _map.Vol, _map.Vol, _map.MaxTauBin, _map.MaxTauBin], dtype=np.complex)

    Wout = np.zeros((Wtot.shape[2],1,Wtot.shape[3],1), dtype=np.complex64)
    Wout[:,0,:,0] = Wtot[UPUP, UPUP, :, :] # r1, r2=0, t1, t2=0

    Win = np.zeros((1,Wtot.shape[2],1,Wtot.shape[3]), dtype=np.complex64)
    Win[0,:,0,:] = Wtot[UPUP, UPUP, :, :] # r1=0, r2, t1=0, t2

    GammaW1=FitGammaW(GammaW, _map)
    GammaW2=RestoreGammaW(GammaW1, _map)
    # if np.allclose(GammaW, GammaW2, rtol=1e-3, atol=1e-5):
        # print "Basis works well!"
    # else:
        # print "Basis fails to work!"
    # print np.amax(np.abs(GammaW-GammaW2))
    # print (GammaW[1,0,0,:,:]).diagonal()
    # print (GammaW2[1,0,0,:,:]).diagonal()
    # print (GammaW[1,0,0,:,:]-GammaW2[1,0,0,:,:]).diagonal()

    print (GammaW[1,0,0,15,:])
    print (GammaW2[1,0,0,15,:])
    print (GammaW[1,0,0,15,:]-GammaW2[1,0,0,15,:])
    for r1 in range(8):
        print "0", r1, np.amax(np.abs(GammaW[0,r1,:,:]-GammaW2[0,r1,:,:]))
    for r1 in range(8):
        print "1", r1,t, np.amax(np.abs(GammaW[1,r1,0,:,:]-GammaW2[1,r1,0,:,:]))
    # for r1 in range(_map.Vol):
        # for t in range(_map.MaxTauBin):
            # print "1", r1,t, np.amax(np.abs(GammaW[1,r1,0,t,:]-GammaW2[1,r1,0,t,:]))
    GammaW=GammaW2

    #Type 0 and 1
    for s in range(2):
        print "Calculate GammaW {0}".format(s)
        TempGammaW=FFTGammaW(GammaW[s,...], _map, 1)
        print "Calculate GammaW FFT done"
        # WWGammaW=Wout*WWGammaW*Win
        if s==0:
            TempGammaW*=Wout
            TempGammaW*=Win
        else:
            #s==1
            TempGammaW*=Wout*2
            TempGammaW*=Win*2

        # TempGammaW=Wout*TempGammaW*Win
        print "Calculate GammaW FFT back"
        WWGammaW[s,...]=FFTGammaW(TempGammaW, _map, -1)
        print "Calculate GammaW done"
        log.info(green("Memory Usage before collecting: {0} MB".format(memory_usage())))
        gc.collect()
        log.info(green("Memory Usage : {0} MB".format(memory_usage())))
    W0.FFT("R","T")
    W.FFT("R","T")

    print "calculating WWGammaW with fourier done!"
    return -1.0*WWGammaW

def AddG_To_WWGammaW(WWGammaW, G, _map):
    #integer tin and tout
    G.FFT("R","T")
    spinUP=_map.Spin2Index(UP,UP)
    spinDOWN=_map.Spin2Index(DOWN,DOWN)
    sub=0
    GammaG = np.zeros([2, _map.Vol, _map.MaxTauBin, _map.MaxTauBin])+0.0*1j

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

                # GammaG[UP, r, tout, tin] += Gout[UP,UP]*WWGammaW[0,r,r,tout,tin]
                # GammaG[UP, r, tout, tin] += Gout[DOWN,DOWN]*WWGammaW[5,r,r,tout,tin]
                # GammaG[DOWN, r, tout, tin] += Gout[DOWN,DOWN]*WWGammaW[1,r,r,tout,tin]
                # GammaG[DOWN, r, tout, tin] += Gout[UP,UP]*WWGammaW[4,r,r,tout,tin]

                GammaG[UP, r, tout, tin] += Gout[UP,UP]*WWGammaW[0,r,r,tout,tin]
                GammaG[UP, r, tout, tin] += Gout[DOWN,DOWN]*WWGammaW[1,r,r,tout,tin]
                GammaG[DOWN, r, tout, tin] += Gout[DOWN,DOWN]*(-np.conj(WWGammaW[0,r,r,tout,tin]))
                GammaG[DOWN, r, tout, tin] += Gout[UP,UP]*(-np.conj(WWGammaW[1,r,r,tout,tin]))
    return GammaG

def shift(r, L):
    if r<0:
        r+=L
    elif r>=L:
        r-=L
    return r

def FullGGGammaG(IrGGGammaG, W0, _map):
    sub=0
    BKPolar=weight.Weight("SmoothT", _map, "FourSpins", "Symmetric","R","T")

    UPUP=_map.Spin2Index(UP,UP)
    DOWNDOWN=_map.Spin2Index(DOWN,DOWN)
    UPDOWN=_map.Spin2Index(UP,DOWN)
    DOWNUP=_map.Spin2Index(DOWN,UP)

    IrGGGammaGuu=np.zeros((_map.Vol, _map.MaxTauBin))+0.0*1j
    IrGGGammaGdu=np.zeros((_map.Vol, _map.MaxTauBin))+0.0*1j
    for i in range(_map.MaxTauBin):
        IrGGGammaGuu[:, i]=IrGGGammaG[UP, :, i, i]
        IrGGGammaGdu[:, i]=IrGGGammaG[DOWN, :, i, i]

    BKPolar.Data[UPUP, sub, UPUP, sub, :,:]=IrGGGammaGuu 
    BKPolar.Data[DOWNDOWN, sub, DOWNDOWN, sub, :,:]=IrGGGammaGuu 
    BKPolar.Data[DOWNDOWN, sub, UPUP, sub, :,:]=IrGGGammaGdu 
    BKPolar.Data[UPUP, sub, DOWNDOWN, sub, :,:]=IrGGGammaGdu 

    # print "BKPolar[UP,UP]=\n", BKPolar.Data[UPUP,0,UPUP,0,0,:]

    BKPolar.FFT("K", "W")
    W0.FFT("K")

    Denorm,JP=calc.Calculate_Denorminator(W0, BKPolar, _map)

    # JPJ=np.einsum("ijklvt,klmnv->ijmnvt", JP, W0.Data)
    BKChiTensor=weight.Weight("SmoothT", _map, "FourSpins", "Symmetric", "K","W")
    lu_piv,Determ=weight.LUFactor(Denorm)
    Check_Denorminator(Denorm, Determ, _map)
    BKChiTensor.LUSolve(lu_piv, -BKPolar.Data)
    return BKChiTensor, Determ


def Calculate_Chi(ChiTensor, _map):
    NSpin, NSublat=ChiTensor.NSpin, ChiTensor.NSublat
    SxSx=np.zeros((NSpin,NSpin))
    SySy=np.zeros((NSpin,NSpin))
    SzSz=np.zeros((NSpin,NSpin))
    UU=_map.Spin2Index(UP,UP) 
    UD=_map.Spin2Index(UP,DOWN) 
    DU=_map.Spin2Index(DOWN, UP) 
    DD=_map.Spin2Index(DOWN, DOWN) 
    SxSx[UD, UD]=SxSx[DU, DU]=1
    SxSx[UD, DU]=SxSx[DU, UD]=1
    SySy[UD, UD]= SySy[DU, DU]=-1
    SySy[UD, DU]= SySy[DU, UD]=1
    SzSz[UU, UU]= SzSz[DD, DD]=1
    SzSz[UU, DD]= SzSz[DD, UU]=-1
    Chi=weight.Weight("SmoothT", _map, "NoSpin", "Symmetric", ChiTensor.SpaceDomain, ChiTensor.TimeDomain)
    Chi.Data=0.0

    SS=[SzSz/4.0]
    for i in range(len(SS)):
        temp=np.einsum("ik, kminvt->mnvt", SS[i], ChiTensor.Data)
        Chi.Data+=temp.reshape([1, NSublat, 1, NSublat, _map.Vol, _map.MaxTauBin]) 
    return Chi

def Check_Denorminator(Denorm, Determ, _map):
    pos=np.where(Determ==Determ.min())
    x,t=pos[0][0], pos[1][0]
    log.info("Checking denorminator of GammaG")
    log.info("The minmum {0} is at K={1} and Omega={2}".format(Determ.min(), _map.IndexToCoordi(x), t))
    SpSub,Vol,Time=Denorm.shape[0]*Denorm.shape[1], Denorm.shape[-2], Denorm.shape[-1]
    Denorm=Denorm.reshape([SpSub,SpSub,Vol,Time])
    log.info("The 1/linalg.cond is {0}".format(1.0/np.linalg.cond(Denorm[...,x,t])))
    if Determ.min().real<0.0 and Determ.min().imag<1.0e-4:
        raise DenorminatorTouchZero(Determ.min(), _map.IndexToCoordi(x), t)

def FitGammaW(GammaW, _map):
    BasisNum=32
    TauBin=_map.MaxTauBin
    svd=basis.SVDBasis(TauBin, _map.Beta, "Bose")
    svd.GenerateBasis(BasisNum)
    Basis=svd.GetBasis()["Basis"]
    # for e in range(BasisNum):
        # Basis[e]=Basis[e][::4]
    BasisVec=np.array(Basis)  #BasisNum, MaxTauBin
    NewGammaW=np.einsum("sijmn,km->sijkn", GammaW, BasisVec) 
    NewGammaW=np.einsum("sijkn,ln->sijkl", NewGammaW, BasisVec) 
    return NewGammaW

def RestoreGammaW(GammaW, _map):
    BasisNum=32
    TauBin=_map.MaxTauBin
    svd=basis.SVDBasis(TauBin, _map.Beta, "Bose")
    svd.GenerateBasis(BasisNum)
    Basis=svd.GetBasis()["Basis"]
    # for e in range(BasisNum):
        # Basis[e]=Basis[e][::4]
    BasisVec=np.array(Basis)  #BasisNum, MaxTauBin
    NewGammaW=np.einsum("sijkl,km->sijml", GammaW, BasisVec) 
    NewGammaW=np.einsum("sijml,ln->sijmn", NewGammaW, BasisVec) 
    return NewGammaW

