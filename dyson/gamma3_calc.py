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
            GGammaG[UP, :, t, :]+=Gout[UP,UP]*GammaG[UP,:,tout,:]
            GGammaG[DOWN, :, t, :]+=Gout[DOWN,DOWN]*GammaG[DOWN,:,tout,:]

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
            GGGammaG[UP, :, :, t]+=Gin[UP,UP]*GGammaG[UP,:,:,tin]
            GGGammaG[DOWN, :, :, t]+=Gin[DOWN,DOWN]*GGammaG[DOWN,:,:,tin]
    return GGGammaG*_map.Beta**2/_map.MaxTauBin**2

def AddG_To_GGGammaG(GGGammaG, G, _map):
    #integer tin and tout
    G.FFT("R","T")
    FermiLoopSign=-1
    spinUP=_map.Spin2Index(UP,UP)
    spinDOWN=_map.Spin2Index(DOWN,DOWN)
    sub=0
    #only GammaW(r,r) has non-zero values, so that we only keep one r for now
    LocalGammaW = np.zeros([2, _map.Vol, _map.MaxTauBin, _map.MaxTauBin], dtype=np.complex)
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

            ## UP UP UP UP
            #GammaW[0, r, r, tout, tin] = Gout[UP,UP]*GGGammaG[UP,r,tout,tin]

            ## DOWN DOWN DOWN DOWN
            #GammaW[1, r, r, tout, tin] = Gout[DOWN,DOWN]*GGGammaG[DOWN,r,tout,tin]

            ## in:UP DOWN out:DOWN UP
            #GammaW[5, r, r, tout, tin] = Gout[UP, UP]*GGGammaG[DOWN,r,tout,tin]

            ## in:DOWN UP out:UP DOWN 
            #GammaW[4, r, r, tout, tin] = Gout[DOWN,DOWN]*GGGammaG[UP,r,tout,tin]

            # UP UP UP UP+ DOWN DOWN DOWN DOWN
            LocalGammaW[0, :, tout, tin] = Gout[UP,UP]*GGGammaG[UP,:,tout,tin]+Gout[DOWN,DOWN]*GGGammaG[DOWN,:,tout,tin]

            # in:DOWN UP out:UP DOWN 
            LocalGammaW[1, :, tout, tin] = Gout[DOWN,DOWN]*GGGammaG[UP,:,tout,tin]

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

            ## UP UP UP UP
            #GammaW[0, r, r, tout, tin] += Gout[UP,UP]*GGGammaG[UP,r,tin,tout]

            ## DOWN DOWN DOWN DOWN
            #GammaW[1, r, r, tout, tin] += Gout[DOWN,DOWN]*GGGammaG[DOWN,r,tin,tout]

            ## in:UP DOWN out:DOWN UP
            #GammaW[5, r, r, tout, tin] += Gout[DOWN, DOWN]*GGGammaG[UP,r,tin,tout]

            ## in:DOWN UP out:UP DOWN 
            #GammaW[4, r, r, tout, tin] += Gout[UP,UP]*GGGammaG[DOWN,r,tin,tout]

            ## UP UP UP UP + DOWN DOWN DOWN DOWN
            LocalGammaW[0, :, tout, tin] += Gout[UP,UP]*GGGammaG[UP,:,tin,tout]+Gout[DOWN,DOWN]*GGGammaG[DOWN,:,tin,tout]

            ## in:DOWN UP out:UP DOWN 
            LocalGammaW[1, :, tout, tin] += Gout[UP,UP]*GGGammaG[DOWN,:,tin,tout]

    GammaW = np.zeros([2, _map.Vol, _map.Vol, _map.MaxTauBin, _map.MaxTauBin], dtype=np.complex)
    for r in range(_map.Vol):
        GammaW[:, r, r, :, :]=LocalGammaW[:,r,:,:]

    # GammaWDiagonal=FitGammaW_diagonal(GammaW, _map)
    # GammaW = np.zeros([2, _map.Vol, _map.Vol, _map.BasisNum, _map.BasisNum], dtype=np.complex)
    # for r in range(_map.Vol):
        # GammaW[:, r, r, :, :]=GammaWDiagonal[:,r,:,:]

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

def FFTGammaW_Space(GammaW, _map, BackForth):
    OldShape=GammaW.shape
    NewShape=(_map.L[0], _map.L[1], _map.L[0], _map.L[1], _map.BasisNum, _map.BasisNum)
    GammaW=GammaW.reshape(NewShape)
    if BackForth==1:
        GammaW=np.fft.fftn(GammaW, axes=[0,1,2,3]) 
    elif BackForth==-1:
        GammaW=np.fft.ifftn(GammaW, axes=[0,1,2,3]) 
    GammaW=GammaW.reshape(OldShape)
    return GammaW

def FFTWshift_Space(Wshift, _map, BackForth):
    OldShape=Wshift.shape
    NewShape=(_map.L[0], _map.L[1], _map.MaxTauBin, _map.MaxTauBin)
    Wshift=Wshift.reshape(NewShape)
    if BackForth==1:
        Wshift=np.fft.fftn(Wshift, axes=[0,1]) 
    elif BackForth==-1:
        Wshift=np.fft.ifftn(Wshift, axes=[0,1]) 
    Wshift=Wshift.reshape(OldShape)
    return Wshift


def AddTwoW_To_GammaW_basis(GammaW, W0, W, _map):
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
    Wshift=np.array(Wshift.Data[:,0,:,0,:,:])*_map.Beta/_map.MaxTauBin
    Wtot=np.zeros([_map.Vol, _map.MaxTauBin, _map.MaxTauBin], dtype=np.complex)
    Wshift[:,:,:,0]+=W0.Data[:,0,:,0,:]
    for t1 in range(_map.MaxTauBin):
        for t2 in range(_map.MaxTauBin):
            dt=t1-t2
            if dt<0:
                dt+=_map.MaxTauBin
            Wtot[:, t1, t2]=Wshift[0,0,:,dt]
    Wtot=FFTWshift_Space(Wtot, _map, 1)
    Wtot=FitW(Wtot, _map)

    GammaW1=RestoreGammaW(GammaW, _map)
    print "GammaW"
    print GammaW1[0,0,0,0,:]
    print GammaW1[1,0,0,15,:]

    WWGammaW=np.zeros([2, _map.Vol, _map.Vol, _map.BasisNum, _map.BasisNum], dtype=np.complex)

    # FittedGammaW=FitGammaW(GammaW, _map)
    for s in range(2):
        TempGammaW=FFTGammaW_Space(GammaW[s,...], _map, 1)
        TempGammaW=np.einsum("ijkl, imk->ijml", TempGammaW, Wtot)
        TempGammaW=np.einsum("ijml, jnl->ijmn", TempGammaW, Wtot)
        TempGammaW=FFTGammaW_Space(TempGammaW, _map, -1)
        if s==1:
            TempGammaW*=4.0
        WWGammaW[s,...]=TempGammaW

    # WWGammaW=RestoreGammaW(WWGammaW, _map)

    log.info(green("Memory Usage before collecting: {0} MB".format(memory_usage())))
    gc.collect()
    log.info(green("Memory Usage : {0} MB".format(memory_usage())))
    W0.FFT("R","T")
    W.FFT("R","T")

    return -1.0*WWGammaW

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

    # W.Data=(W.Data+W.Data[:,:,:,:,:,::-1])/2
    #when necessary, one may need to symmetrize W

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

    Wout = np.zeros((Wtot.shape[2],1,Wtot.shape[3],1), dtype=np.complex)
    Wout[:,0,:,0] = Wtot[UPUP, UPUP, :, :] # r1, r2=0, t1, t2=0

    Win = np.zeros((1,Wtot.shape[2],1,Wtot.shape[3]), dtype=np.complex)
    Win[0,:,0,:] = Wtot[UPUP, UPUP, :, :] # r1=0, r2, t1=0, t2

    GammaW=UnCompressGammaW(GammaW, _map)

    # print "GammaW"
    # # GammaW1=FitGammaW(GammaW, _map)
    # # GammaW2=RestoreGammaW(GammaW1, _map)
    # GammaW1=CompressGammaW(GammaW, _map)
    # GammaW2=UnCompressGammaW(GammaW1, _map)
    # # if np.allclose(GammaW, GammaW2, rtol=1e-3, atol=1e-5):
        # # print "Basis works well!"
    # # else:
        # # print "Basis fails to work!"
    # # print np.amax(np.abs(GammaW-GammaW2))
    # # print (GammaW[1,0,0,:,:]).diagonal()
    # # print (GammaW2[1,0,0,:,:]).diagonal()
    # # print (GammaW[1,0,0,:,:]-GammaW2[1,0,0,:,:]).diagonal()

    # print "1,7"
    # print (GammaW[0,0,0,0,:])
    # print (GammaW2[0,0,0,0,:])
    # # print (GammaW[0,1,7,0,:])
    # # print (GammaW2[0,1,7,0,:])
    # print (GammaW[0,0,0,0,:]-GammaW2[0,0,0,0,:])

    # for r1 in range(8):
        # # for r2 in range(8):
        # print "0", r1, np.amax(np.abs(GammaW[0,r1,:,:,:]-GammaW2[0,r1,:,:,:]))
    # for r1 in range(8):
        # print "1", r1,t, np.amax(np.abs(GammaW[1,r1,:,:,:]-GammaW2[1,r1,:,:,:]))
    # for r1 in range(_map.Vol):
        # for t in range(_map.MaxTauBin):
            # print "1", r1,t, np.amax(np.abs(GammaW[1,r1,0,t,:]-GammaW2[1,r1,0,t,:]))

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

    print "WWGammaW Compress"
    WWGammaW1=CompressGammaW(WWGammaW, _map)
    WWGammaW2=UnCompressGammaW(WWGammaW1, _map)
    for r1 in range(8):
        print "0", r1, np.amax(np.abs(WWGammaW[0,r1,:,:]-WWGammaW2[0,r1,:,:]))
    for r1 in range(8):
        print "1", r1,t, np.amax(np.abs(WWGammaW[1,r1,:,:,:]-WWGammaW2[1,r1,:,:,:]))
    WWGammaW=WWGammaW2
    # # WWGammaW1=FitGammaW(WWGammaW, _map)
    # # WWGammaW2=RestoreGammaW(WWGammaW1, _map)
    WWGammaW1=CompressGammaW1(WWGammaW, _map)
    print WWGammaW1.shape
    print "WWGammaW UnCompress"
    WWGammaW2=UnCompressGammaW1(WWGammaW1, _map)
    print "WWGammaW restored"
    print WWGammaW2.shape
    # # if np.allclose(GammaW, GammaW2, rtol=1e-3, atol=1e-5):
        # # print "Basis works well!"
    # # else:
        # # print "Basis fails to work!"
    # # print np.amax(np.abs(GammaW-GammaW2))
    # # print (GammaW[1,0,0,:,:]).diagonal()
    # # print (GammaW2[1,0,0,:,:]).diagonal()
    # # print (GammaW[1,0,0,:,:]-GammaW2[1,0,0,:,:]).diagonal()

    print WWGammaW[1,1,0,3,8]
    print WWGammaW[1,1,0,32-3,32-8]

    print "WWGammaW"
    print (WWGammaW[1,0,1,2,:])
    # print (WWGammaW[1,0,1,:,1])
    # print (WWGammaW[1,1,0,15,15])
    # print "WWGammaW2"
    print (WWGammaW2[1,0,1,2,:])
    # print (WWGammaW2[1,0,1,:,1])

    print (WWGammaW[1,0,1,2,:]-WWGammaW2[1,0,1,2,:])
    for r1 in range(8):
        print "0", r1, np.amax(np.abs(WWGammaW[0,r1,:,:]-WWGammaW2[0,r1,:,:]))
    # for r1 in range(8):
        # for r2 in range(8):
    r1,r2=0,1
    for t in range(_map.MaxTauBin):
        print "1", r1, r2, t, np.amax(np.abs(WWGammaW[1,:,:,t,:]-WWGammaW2[1,:,:,t,:]))

    return -1.0*WWGammaW

def AddG_To_WWGammaW(WWGammaW, G, _map):
    #integer tin and tout
    G.FFT("R","T")
    spinUP=_map.Spin2Index(UP,UP)
    spinDOWN=_map.Spin2Index(DOWN,DOWN)
    sub=0
    GammaG = np.zeros([2, _map.Vol, _map.MaxTauBin, _map.MaxTauBin])+0.0*1j
    #only GammaW(r,r) is needed!

    WWGammaW=np.einsum("siikl->sikl", WWGammaW)
    # WWGammaW=RestoreGammaW_diagonal(WWGammaW, _map)

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

            GammaG[UP, :, tout, tin] += Gout[UP,UP]*WWGammaW[0,:,tout,tin]
            GammaG[UP, :, tout, tin] += Gout[DOWN,DOWN]*WWGammaW[1,:,tout,tin]
            GammaG[DOWN, :, tout, tin] += Gout[DOWN,DOWN]*(-np.conj(WWGammaW[0,:,tout,tin]))
            GammaG[DOWN, :, tout, tin] += Gout[UP,UP]*(-np.conj(WWGammaW[1,:,tout,tin]))
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

def GetBoseBasis(_map, BasisNum=None):
    TauBin=_map.MaxTauBin
    svd=basis.SVDBasis(TauBin, _map.Beta, "Bose")
    if BasisNum==None:
        BasisNum=_map.BasisNum
    svd.GenerateBasis(BasisNum)
    Basis=svd.GetBasis()["Basis"]
    # for e in range(BasisNum):
        # Basis[e]=Basis[e][::4]
    BasisVec=np.array(Basis)  #BasisNum, MaxTauBin
    return BasisVec

def FitW(Wshift, _map):
    BasisVec=GetBoseBasis(_map)
    NewWshift=np.einsum("imn,km->ikn", Wshift, BasisVec) 
    NewWshift=np.einsum("ikn,ln->ikl", NewWshift, BasisVec) 
    return NewWshift

def FitGammaW_diagonal(GammaW, _map):
    BasisVec=GetBoseBasis(_map)
    NewGammaW=np.einsum("srmn,km->srkn", GammaW, BasisVec) 
    NewGammaW=np.einsum("srkn,ln->srkl", NewGammaW, BasisVec) 
    return NewGammaW

def FitGammaW(GammaW, _map):
    BasisVec=GetBoseBasis(_map)
    NewGammaW=np.einsum("sijmn,km->sijkn", GammaW, BasisVec) 
    NewGammaW=np.einsum("sijkn,ln->sijkl", NewGammaW, BasisVec) 
    return NewGammaW

def RestoreGammaW_diagonal(GammaW, _map):
    BasisVec=GetBoseBasis(_map)
    NewGammaW=np.einsum("srkl,km->srml", GammaW, BasisVec) 
    NewGammaW=np.einsum("srml,ln->srmn", NewGammaW, BasisVec) 
    return NewGammaW

def RestoreGammaW(GammaW, _map, TauBin=None):
    if TauBin==None:
        TauBin=_map.MaxTauBin
    Interval=int(_map.MaxTauBin/TauBin)
    BasisVec=GetBoseBasis(_map)
    ReducedBasisVec=np.zeros([_map.BasisNum, TauBin], dtype=BasisVec.dtype)
    for e in range(_map.BasisNum):
        ReducedBasisVec[e]=BasisVec[e][::Interval]
    NewGammaW=np.einsum("sijkl,km->sijml", GammaW, ReducedBasisVec) 
    NewGammaW=np.einsum("sijml,ln->sijmn", NewGammaW, ReducedBasisVec) 
    return NewGammaW

def CompressGammaW(GammaW, _map):
    if len(GammaW.shape)==2:
        #do not need to compress
        return GammaW
    N=_map.Vol*_map.MaxTauBin
    CompactGammaW=np.zeros((2, N*(N+1)/2) , dtype=complex)
    GammaW=GammaW.swapaxes(2,3)
    GammaW=GammaW.reshape([2, N, N])
    GammaW[0,:]+=-np.conj(GammaW[0,:].T)
    GammaW[1,:]+=-np.conj(GammaW[1,:].T)
    GammaW/=2.0
    Index=np.triu_indices(N)
    CompactGammaW[0,:]=GammaW[0,:][Index]
    CompactGammaW[1,:]=GammaW[1,:][Index]
    return CompactGammaW

def UnCompressGammaW(CompactGammaW, _map):
    if len(CompactGammaW.shape)==5:
        #do not need to uncompress
        return CompactGammaW
    N=_map.Vol*_map.MaxTauBin
    GammaW=np.zeros((2, N, N), dtype=complex)
    uIndex=np.triu_indices(N) #-1 is used to avoid the diagonal elements
    GammaW[0,:][uIndex]=CompactGammaW[0, :]
    GammaW[1,:][uIndex]=CompactGammaW[1, :]

    diag=np.array(np.diag(GammaW[0,:]))
    GammaW[0,:]+=-np.conj(GammaW[0,:].T)
    np.fill_diagonal(GammaW[0,:], diag)
    diag=np.array(np.diag(GammaW[1,:]))
    GammaW[1,:]+=-np.conj(GammaW[1,:].T)
    np.fill_diagonal(GammaW[1,:], diag)

    GammaW=GammaW.reshape([2, _map.Vol, _map.MaxTauBin, _map.Vol, _map.MaxTauBin])
    GammaW=GammaW.swapaxes(2,3)
    return GammaW

# def SymmetryMapping(_map, LatName):
    # # if LatName=="Triangular":
    # Lx, Ly=_map.L[0], _map.L[1]
    # Vol=_map.Vol
    # TauBin=_map.MaxTauBin
    # Squeeze=np.arange(Vol**2*TauBin**2, dtype=int).reshape([Vol, Vol, TauBin, TauBin])
    # Restore=[]
    # Flag=np.zeros([Vol, Vol, TauBin, TauBin], dtype=bool)
    # # Restore=np.zeros([Vol, Vol, TauBin, TauBin], dtype=int)
    # Index=0
    # for r1 in range(Vol):
        # for r2 in range(Vol):
            # for t1 in range(TauBin):
                # for t2 in range(t1, TauBin):
                    # Squeeze[r1,r2,t1,t2]=Index
                    # Restore.append(r1*Vol*TauBin**2+r2*TauBin**2+t1*TauBin+t2)
                    # if t1 is not t2:
                        # Squeeze[r1,r2,t2,t1]=Index
                        # Flag[r1,r2,t2,t1]=True
                    # Index+=1
    # print "Index", Index
    # return Squeeze, Restore, Flag

def GetR(r, L):
    if r<0:
        r+=L
    if r<0:
        r+=L
    if r>=L:
        r-=L
    if r>=L:
        r-=L
    return r

def MirrorR(vec, _map, LatName):
    x1, y1, x2, y2=vec
    Lx,Ly=_map.L[0], _map.L[1]
    Rlist=[vec,]
    
    if LatName=="Triangular":
        Rlist.append((GetR(-x1-y1, Lx), y1, GetR(-x2-y2, Lx), y2)) 
        Rlist.append((GetR(-x1, Lx), GetR(-y1, Ly), GetR(-x2, Lx), GetR(-y2, Ly))) 
        Rlist.append((GetR(x1+y1, Lx), GetR(-y1, Ly), GetR(x2+y2, Lx), GetR(-y2, Ly))) 
    elif LatName=="Square":
        Rlist.append((GetR(-x1, Lx), y1, GetR(-x2, Lx), y2)) 
        Rlist.append((GetR(-x1, Lx), GetR(-y1, Ly), GetR(-x2, Lx), GetR(-y2, Ly))) 
        Rlist.append((GetR(x1, Lx), GetR(-y1, Ly), GetR(x2, Lx), GetR(-y2, Ly))) 
    return set(Rlist)

def SymmetryMapping(_map, LatName):
    # if LatName=="Triangular":
    print "test:", MirrorR((1,0,0,1), _map, "Triangular")
    print "test:", MirrorR((1,0,7,2), _map, "Triangular")
    print "test:", MirrorR((1,0,1,1), _map, "Triangular")
    Lx, Ly=_map.L[0], _map.L[1]
    Vol=_map.Vol
    TauBin=_map.MaxTauBin
    TauSqueeze=np.zeros([TauBin, TauBin], dtype=int)
    TauFlag=np.zeros([TauBin, TauBin], dtype=int)
    Index=0
    TauRestore=[]
    for t1 in range(TauBin):
        for t2 in range(t1, TauBin):
            TauSqueeze[t1,t2]=Index
            if t1 is not t2:
                TauSqueeze[t2,t1]=Index
                TauFlag[t2,t1]=True
            TauRestore.append((t1, t2))
            Index+=1

    RSqueeze=np.zeros([Vol, Vol], dtype=int)
    RRestore=[]
    # for r1 in range(Vol):
        # for r2 in range(Vol):
    Points=[]
    Index=0
    for x1 in range(Lx):
        for y1 in range(Ly):
            for x2 in range(Lx):
                for y2 in range(Ly):
                    vec=(x1,y1,x2,y2)
                    if vec in Points:
                        continue
                    Rlist=MirrorR(vec, _map, LatName)
                    # Rlist=MirrorR(vec, _map, "Square")
                    # for e in Rlist:
                        # if e in Points:
                            # print "error:", e
                    # print len(Rlist)
                    for e in Rlist:
                        Points.append(e)
                        r1=_map.CoordiIndex((e[0],e[1]))
                        r2=_map.CoordiIndex((e[2],e[3]))
                        RSqueeze[r1,r2]=Index
                    Index+=1
                    e=list(Rlist)[0]
                    r1=_map.CoordiIndex((e[0],e[1]))
                    r2=_map.CoordiIndex((e[2],e[3]))
                    RRestore.append((r1,r2))
    RSize=Index
    return TauSqueeze, TauRestore, TauFlag, RSqueeze, RRestore


def CompressGammaW1(GammaW, _map):
    if len(GammaW.shape)==2:
        #do not need to compress
        return GammaW
    Lx, Ly=_map.L[0], _map.L[1]
    Vol=_map.Vol
    TauBin=_map.MaxTauBin
    TauSqueeze, TauRestore, TauFlag, RSqueeze, RRestore=SymmetryMapping(_map, "Triangular")
    CompactGammaW=np.zeros([2, _map.Vol, _map.Vol, len(TauRestore)], dtype=complex)
    # for s in range(2):
    for t1 in range(TauBin):
        for t2 in range(t1, TauBin):
            CompactGammaW[:,:,:,TauSqueeze[t1,t2]]=GammaW[:,:,:,t1,t2]

    RCompactGammaW=np.zeros([2, len(RRestore), len(TauRestore)], dtype=complex)
    # print RRestore
    for i in range(len(RRestore)):
        r1,r2=RRestore[i]
        # print i, r1, r2
        RCompactGammaW[:,i,:]=CompactGammaW[:,r1,r2,:]
    return RCompactGammaW

def UnCompressGammaW1(CompactGammaW, _map):
    import numpy.ma as npma
    if len(CompactGammaW.shape)==5:
        #do not need to uncompress
        return CompactGammaW
    Lx, Ly=_map.L[0], _map.L[1]
    Vol=_map.Vol
    TauBin=_map.MaxTauBin
    TauSqueeze, TauRestore, TauFlag, RSqueeze, RRestore=SymmetryMapping(_map, "Triangular")

    TauGammaW=np.zeros((2, Vol, Vol, len(TauRestore)), dtype=complex)
    for r1 in range(Vol):
        for r2 in range(Vol):
            TauGammaW[:,r1,r2,:]=CompactGammaW[:,RSqueeze[r1,r2],:]

    GammaW=np.zeros((2, Vol, Vol, TauBin, TauBin), dtype=complex)
    for t1 in range(TauBin):
        for t2 in range(0, t1):
            GammaW[:,:,:,t1,t2]=-np.conj(TauGammaW[:,:,:,TauSqueeze[t1,t2]])
    GammaW=GammaW.swapaxes(1,2)
    for t1 in range(TauBin):
        for t2 in range(t1, TauBin):
            GammaW[:,:,:,t1,t2]=TauGammaW[:,:,:,TauSqueeze[t1,t2]]
    return GammaW


# def CompressGammaW(GammaW, _map):
    # Interval=int(_map.MaxTauBin/_map.MaxTauBinTiny)
    # GammaW=GammaW[:,:,:,::Interval,::Interval]
    # return GammaW

# def UnCompressLocalGammaW(GammaW, _map):
    # BasisNum=32
    # GammaW=np.einsum("siikl->sikl", GammaW)
    # BasisVec=GetBoseBasis(_map, BasisNum)
    # Interval=int(_map.MaxTauBin/_map.MaxTauBinTiny)
    # ReducedBasisVec=np.zeros([BasisNum, _map.MaxTauBinTiny], dtype=BasisVec.dtype)
    # for e in range(BasisNum):
        # ReducedBasisVec[e]=BasisVec[e][::Interval]
    # NewGammaW=np.einsum("srmn,km->srkn", GammaW, ReducedBasisVec) 
    # NewGammaW=np.einsum("srkn,ln->srkl", NewGammaW, ReducedBasisVec) 
    # NewGammaW*=Interval**2

    # NewGammaW=np.einsum("sikl,km->siml", NewGammaW, BasisVec) 
    # NewGammaW=np.einsum("siml,ln->simn", NewGammaW, BasisVec) 
    # return NewGammaW

# def MirrorMapping(_map):
    # Lx, Ly=_map.L[0], _map.L[1]
    # RestoreMap=np.zeros((Lx/2+1)*(Ly/2+1), dtype=int)
    # MirrorMap=np.zeros(Lx*Ly, dtype=int)
    # for x in range(Lx):
        # for y in range(Ly):
            # xp,yp=x,y
            # if x>Lx/2:
                # xp=Lx-xp
            # if y>Ly/2:
                # yp=Ly-yp
            # MirrorMap[_map.CoordiIndex([x,y])]=xp*(Ly/2+1)+yp

    # for xp in range(Lx/2+1):
        # for yp in range(Ly/2+1):
            # x,y=xp,yp
            # RestoreMap[xp*(Ly/2+1)+yp]=_map.CoordiIndex([x,y])
    # print MirrorMap, RestoreMap

    # return (MirrorMap, RestoreMap)

# def CompressGammaW_MirrorSymmetry(GammaW, _map):
    # Lx, Ly=_map.L[0], _map.L[1]
    # SmallVol=(Lx/2+1)*(Ly/2+1)
    # TauBin=GammaW.shape[-1]
    # NewGammaW=np.zeros([2, SmallVol, SmallVol, TauBin, TauBin], dtype=np.complex)
    # MirrorMap, RestoreMap=MirrorMapping(_map)
    # for r1p in range(SmallVol):
        # for r2p in range(SmallVol):
            # r1=RestoreMap[r1p]
            # r2=RestoreMap[r2p]
            # # print (r1p,r2p, r1, r2)
            # NewGammaW[:,r1p,r2p,:,:]=GammaW[:,r1,r2,:,:]
    # return NewGammaW

# def UnCompressGammaW_MirrorSymmetry(GammaW, _map):
    # Lx, Ly=_map.L[0], _map.L[1]
    # SmallVol=(Lx/2+1)*(Ly/2+1)
    # TauBin=GammaW.shape[-1]
    # NewGammaW=np.zeros([2, _map.Vol, _map.Vol, TauBin, TauBin], dtype=np.complex)
    # MirrorMap, RestoreMap=MirrorMapping(_map)
    # for r1 in range(_map.Vol):
        # for r2 in range(_map.Vol):
            # r1p=MirrorMap[r1]
            # r2p=MirrorMap[r2]
            # NewGammaW[:,r1,r2,:,:]=GammaW[:,r1p,r2p,:,:]
    # return NewGammaW




