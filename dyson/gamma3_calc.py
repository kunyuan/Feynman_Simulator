#!usr/bin/env python
import sys
import numpy as np
import parameter as para
from weight import UP,DOWN,IN,OUT
import weight, plot
from logger import *
import r_index 

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
            # if t1==t2:
                # print W.Data[spinUP,sub,spinUP,sub,0,abs(t1-t2)]
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
            Gout=0.5*signout*G.Data[:,sub,:,sub,r,tgout]
            tgout=t-tout
            signout=1
            if tgout<0:
                tgout+=_map.MaxTauBin
                signout*=-1
            Gout+=0.5*signout*G.Data[:,sub,:,sub,r,tgout]
            for tin in range(_map.MaxTauBin):
                GGGammaG[UP, :, t, tin]+=Gout[UP,UP]*GammaG[UP,:,tout,tin]
                GGGammaG[DOWN, :, t, tin]+=Gout[DOWN,DOWN]*GammaG[DOWN,:,tout,tin]

    GGGammaGNew=np.zeros([2, _map.Vol, _map.MaxTauBin, _map.MaxTauBin])+0.0*1j
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
                GGGammaGNew[UP, :, tout, t]+=Gin[UP,UP]*GGGammaG[UP,:,tout,tin]
                GGGammaGNew[DOWN, :, tout, t]+=Gin[DOWN,DOWN]*GGGammaG[DOWN,:,tout,tin]
    return GGGammaGNew*_map.Beta**2/_map.MaxTauBin**2

def AddG_To_GammaG(GammaG, G, _map):
    #integer tin and tout
    G.FFT("R","T")
    FermiLoopSign=-1
    spinUP=_map.Spin2Index(UP,UP)
    spinDOWN=_map.Spin2Index(DOWN,DOWN)
    sub=0
    GGammaG = np.zeros([6, _map.Vol, _map.Vol, _map.MaxTauBin, _map.MaxTauBin])+0.0*1j
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

                ## UP UP UP UP
                GGammaG[0, r, r, tout, tin] = Gout[UP,UP]*GammaG[UP,r,tout,tin]

                ## DOWN DOWN DOWN DOWN
                GGammaG[1, r, r, tout, tin] = Gout[DOWN,DOWN]*GammaG[DOWN,r,tout,tin]

                ## in:UP DOWN out:DOWN UP
                GGammaG[5, r, r, tout, tin] = Gout[UP, UP]*GammaG[DOWN,r,tout,tin]

                ## in:DOWN UP out:UP DOWN 
                GGammaG[4, r, r, tout, tin] = Gout[DOWN,DOWN]*GammaG[UP,r,tout,tin]

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
                GGammaG[0, r, r, tout, tin] += Gout[UP,UP]*GammaG[UP,r,tin,tout]

                ## DOWN DOWN DOWN DOWN
                GGammaG[1, r, r, tout, tin] += Gout[DOWN,DOWN]*GammaG[DOWN,r,tin,tout]

                ## in:UP DOWN out:DOWN UP
                GGammaG[5, r, r, tout, tin] += Gout[DOWN, DOWN]*GammaG[UP,r,tin,tout]

                ## in:DOWN UP out:UP DOWN 
                GGammaG[4, r, r, tout, tin] += Gout[UP,UP]*GammaG[DOWN,r,tin,tout]
    return FermiLoopSign*GGammaG

def GenerateSpinIndex(_map):
    UPUP=_map.Spin2Index(UP,UP)
    DOWNDOWN=_map.Spin2Index(DOWN,DOWN)
    UPDOWN=_map.Spin2Index(UP,DOWN)
    DOWNUP=_map.Spin2Index(DOWN,UP)

    spinindex=np.array([UP,DOWN])
    spin2index=np.array([UPUP,DOWNDOWN, UPDOWN, DOWNUP])
    return (spinindex, spin2index)

def FastWWGammaW(GGammaG, W0, W, _map):
    import gamma3
    sub = 0
    W.FFT("R","T")
    W0.FFT("R")

    rIndex=np.zeros([_map.Vol, _map.Vol])
    for r in range(_map.Vol):
        for rout in range(_map.Vol):
            rIndex[r, rout] = int(r_index.CoordiIndex(r, rout, _map))

    spinindex, spin2index=GenerateSpinIndex(_map)
    # print "before", GGammaG[:,0,0,0,0]
    WWGammaW=gamma3.fast_wwgammaw(GGammaG, W0.Data[:,sub,:,sub,:], W.Data[:,sub,:,sub,:,:], _map.Beta, rIndex, spinindex, spin2index, _map.Vol, _map.MaxTauBin)
    # print WWGammaW[:,0,0,0,0]
    # print "after", GGammaG[:,0,0,0,0]
    return WWGammaW

def FFTGammaW(GammaW, _map, BackForth):
    OldShape=GammaW.shape
    NewShape=(6, _map.L[0], _map.L[1], _map.L[0], _map.L[1], _map.MaxTauBin, _map.MaxTauBin)
    GammaW=GammaW.reshape(NewShape)
    if BackForth==1:
        GammaW=np.fft.fftn(GammaW, axes=[1,2,3,4,5,6]) 
    elif BackForth==-1:
        GammaW=np.fft.ifftn(GammaW, axes=[1,2,3,4,5,6]) 
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

def FourierWWGammaW(GGammaG, W0, W, _map):
    # import gamma3
    GGammaG=np.array(GGammaG)
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

    GGammaG=FFTGammaW(GGammaG, _map, 1)

    WGammaW=np.zeros([6, _map.Vol, _map.Vol, _map.MaxTauBin, _map.MaxTauBin])+0.0*1j
    Wout=Wtot
    print "calculating WGammaW with fourier..."
    for kout in range(_map.Vol):
            for wout in range(_map.MaxTauBin):
                    # UPUP UPUP
                    WGammaW[0, kout, :, wout, :]  = Wout[UPUP, UPUP, kout, wout] * GGammaG[0, kout, :, wout, :]

                    # DOWNDOWN DOWNDOWN
                    WGammaW[1, kout, :, wout, :]  = Wout[DOWNDOWN, DOWNDOWN, kout, wout] * GGammaG[1, kout, :, wout, :]

                    # out:UPUP in:DOWNDOWN 
                    WGammaW[2, kout, :, wout, :]  = Wout[UPUP, DOWNDOWN, kout, wout] * GGammaG[1, kout, :, wout, :]

                    # out:DOWNDOWN in:UPUP
                    WGammaW[3, kout, :, wout, :]  = Wout[DOWNDOWN, UPUP, kout, wout] * GGammaG[0, kout, :, wout, :]

                    # out:UPDOWN in:DOWNUP
                    WGammaW[4, kout, :, wout, :]  = Wout[UPDOWN, DOWNUP, kout, wout] * GGammaG[5, kout, :, wout, :]

                    # out:DOWNUP in:UPDOWN
                    WGammaW[5, kout, :, wout, :]  = Wout[DOWNUP, UPDOWN, kout, wout] * GGammaG[4, kout, :, wout, :]

    Win=Wtot
    print "calculating WWGammaW with fourier..."

    WWGammaW = np.zeros([6, _map.Vol, _map.Vol, _map.MaxTauBin, _map.MaxTauBin]) + 0.0*1j
    for kin in range(_map.Vol):
            for win in range(_map.MaxTauBin):
                    ## out:UPUP in:UPUP
                    WWGammaW[0, :, kin, :, win] = Win[UPUP, UPUP, kin, win] * WGammaW[0, :, kin, :, win] + Win[UPUP, DOWNDOWN,kin, win] * WGammaW[2, :, kin, :, win]

                    ## out:DOWNDOWN in:DOWNDOWN
                    WWGammaW[1, :, kin, :, win] = Win[DOWNDOWN, UPUP,kin, win] * WGammaW[3, :, kin, :, win] + Win[DOWNDOWN, DOWNDOWN,kin, win] * WGammaW[1, :, kin, :, win]

                    ## out:UPUP in:DOWNDOWN
                    WWGammaW[2, :, kin, :, win] = Win[DOWNDOWN, UPUP,kin, win] * WGammaW[0, :, kin, :, win] + Win[DOWNDOWN, DOWNDOWN,kin, win] * WGammaW[2, :, kin, :, win]

                    ## out:DOWNDOWN in:UPUP
                    WWGammaW[3, :, kin, :, win] = Win[UPUP, UPUP,kin, win] * WGammaW[3, :, kin, :, win] + Win[UPUP, DOWNDOWN, kin, win] * WGammaW[1, :, kin, :, win]

                    ## out:UPDOWN in:DOWNUP
                    WWGammaW[4, :, kin, :, win] = Win[UPDOWN, DOWNUP, kin, win] * WGammaW[4, :, kin, :, win]

                    ## out:DOWNUP in:UPDOWN
                    WWGammaW[5, :, kin, :, win] = Win[DOWNUP, UPDOWN, kin, win] * WGammaW[5, :, kin, :, win]

    WWGammaW=FFTGammaW(WWGammaW, _map, -1)
    W0.FFT("R","T")
    W.FFT("R","T")
    return -1.0*WWGammaW

def FastFourierWWGammaW(GGammaG, W0, W, _map):
    import gamma3
    # GGammaG=np.array(GGammaG)
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

    GGammaG=FFTGammaW(GGammaG, _map, 1)

    WWGammaW=gamma3.fast_fourier_wwgammaw(GGammaG, Wtot, _map.Beta, spinindex, spin2index, _map.Vol, _map.MaxTauBin)

    WWGammaW=FFTGammaW(WWGammaW, _map, -1)
    W0.FFT("R","T")
    W.FFT("R","T")
    return WWGammaW

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
                    Wout = 0.5*W.Data[:,sub,:,sub,dr_out,dt_out]*(_map.Beta/_map.MaxTauBin)
                    dt_out = t - tout
                    if dt_out<0:
                        dt_out+=_map.MaxTauBin
                    Wout += 0.5*W.Data[:,sub,:,sub,dr_out,dt_out]*(_map.Beta/_map.MaxTauBin)

                    if t == tout:
                        Wout += W0.Data[:,sub,:,sub,dr_out]

                    #For Delta W test
                    # if t!= tout:
                        # continue
                    # else:
                        # Wout = W0.Data[:,sub,:,sub,dr_out]

                    # UPUP UPUP
                    WGammaW[0, r, :, t, :]  += Wout[UPUP, UPUP] * GGammaG[0, rout, :, tout, :]

                    # DOWNDOWN DOWNDOWN
                    WGammaW[1, r, :, t, :]  += Wout[DOWNDOWN, DOWNDOWN] * GGammaG[1, rout, :, tout, :]

                    # out:UPUP in:DOWNDOWN 
                    WGammaW[2, r, :, t, :]  += Wout[UPUP, DOWNDOWN] * GGammaG[1, rout, :, tout, :]

                    # out:DOWNDOWN in:UPUP
                    WGammaW[3, r, :, t, :]  += Wout[DOWNDOWN, UPUP] * GGammaG[0, rout, :, tout, :]

                    # out:UPDOWN in:DOWNUP
                    WGammaW[4, r, :, t, :]  += Wout[UPDOWN, DOWNUP] * GGammaG[5, rout, :, tout, :]

                    # out:DOWNUP in:UPDOWN
                    WGammaW[5, r, :, t, :]  += Wout[DOWNUP, UPDOWN] * GGammaG[4, rout, :, tout, :]

    # print "ConventionalWGammaW, type0, tau1=0, dyson=\n", WGammaW[0, 1, 1, 0, :]
    # print "ConventionalWGammaW, type0, tau1=1, dyson=\n", WGammaW[0, 1, 1, 1, :]
    # print "ConventionalWGammaW, type0, tau2=0, dyson=\n", WGammaW[0, 1, 1, :, 0]
    # print "ConventionalWGammaW, type0, tau2=1, dyson=\n", WGammaW[0, 1, 1, :, 1]
    # print "ConventionalWGammaW, type0, diagonal, dyson=\n", WGammaW[0, 1, 1, :, :].diagonal()

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
                    Win = 0.5*W.Data[:,sub,:,sub,dr_in,dt_in]*(_map.Beta/_map.MaxTauBin)
                    dt_in = tin - t-1
                    if dt_in < 0:
                        dt_in +=_map.MaxTauBin
                    Win += 0.5*W.Data[:,sub,:,sub,dr_in,dt_in]*(_map.Beta/_map.MaxTauBin)

                    if t == tin:
                        Win += W0.Data[:,sub,:,sub, dr_in]

                    ##For DeltaW test
                    # if t!= tin:
                        # continue
                    # else:
                        # Win = W0.Data[:,sub,:,sub,dr_in]

                    ## out:UPUP in:UPUP
                    WWGammaW[0, :, r, :, t] += Win[UPUP, UPUP] * WGammaW[0, :, rin, :, tin]
                    WWGammaW[0, :, r, :, t] += Win[UPUP, DOWNDOWN] * WGammaW[2, :, rin, :, tin]

                    ## out:DOWNDOWN in:DOWNDOWN
                    WWGammaW[1, :, r, :, t] += Win[DOWNDOWN, UPUP] * WGammaW[3, :, rin, :, tin]
                    WWGammaW[1, :, r, :, t] += Win[DOWNDOWN, DOWNDOWN] * WGammaW[1, :, rin, :, tin]

                    ## out:UPUP in:DOWNDOWN
                    WWGammaW[2, :, r, :, t] += Win[DOWNDOWN, UPUP] * WGammaW[0, :, rin, :, tin]
                    WWGammaW[2, :, r, :, t] += Win[DOWNDOWN, DOWNDOWN] * WGammaW[2, :, rin, :, tin]

                    ## out:DOWNDOWN in:UPUP
                    WWGammaW[3, :, r, :, t] += Win[UPUP, UPUP] * WGammaW[3, :, rin, :, tin]
                    WWGammaW[3, :, r, :, t] += Win[UPUP, DOWNDOWN] * WGammaW[1, :, rin, :, tin]

                    ## out:UPDOWN in:DOWNUP
                    WWGammaW[4, :, r, :, t] += Win[UPDOWN, DOWNUP] * WGammaW[4, :, rin, :, tin]

                    ## out:DOWNUP in:UPDOWN
                    WWGammaW[5, :, r, :, t] += Win[DOWNUP, UPDOWN] * WGammaW[5, :, rin, :, tin]
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

                GGammaW[UP, r, tout, tin] += Gout[UP,UP]*GammaW[0,r,r,tout,tin]
                GGammaW[UP, r, tout, tin] += Gout[DOWN,DOWN]*GammaW[5,r,r,tout,tin]
                GGammaW[DOWN, r, tout, tin] += Gout[DOWN,DOWN]*GammaW[1,r,r,tout,tin]
                GGammaW[DOWN, r, tout, tin] += Gout[UP,UP]*GammaW[4,r,r,tout,tin]

    GammaG=AddTwoGToGammaG(GGammaW, G, _map)
    return GammaG

def shift(r, L):
    if r<0:
        r+=L
    elif r>=L:
        r-=L
    return r

def FastGammaG_RPA(GammaG, G, W0, _map):
    import gamma3
    #OrderSign=-1, FermiLoopSign=-1, therefore TotalSign=1
    #integer tin and tout
    W0.FFT("R")
    G.FFT("R")
    spinindex, spin2index=GenerateSpinIndex(_map)
    rIndex=np.zeros([_map.Vol, _map.Vol])
    sub=0
    GammaGNew=np.zeros([2, _map.Vol, _map.MaxTauBin, _map.MaxTauBin])+0.0*1j
    for r in range(_map.Vol):
        for rout in range(_map.Vol):
            rIndex[r, rout] = int(r_index.CoordiIndex(r, rout, _map))
    GammaGNew=gamma3.fast_gammag_rpa(GammaGNew, GammaG, G.Data[:,sub,:,sub,0,:], W0.Data[:,sub,:,sub,:], _map.Beta, rIndex, spinindex, spin2index, _map.Vol, _map.MaxTauBin)
    # print GammaGNew
    return GammaGNew


def GammaG_FirstOrder(GammaG, G, W0, _map):
    #OrderSign=-1, FermiLoopSign=-1, therefore TotalSign=1
    #integer tin and tout
    GG=SimpleGG(G,_map)
    GammaGNew=np.zeros([2, _map.Vol, _map.MaxTauBin, _map.MaxTauBin])+0.0*1j
    spinUP=_map.Spin2Index(UP,UP)
    spinDOWN=_map.Spin2Index(DOWN,DOWN)

    # # # print "GammaG[UP,UP]=\n", GammaG[UP,0,:,-1]
    W0.FFT("R")
    G.FFT("R")
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
            # G1=0.5*sign*G.Data[UP,sub,UP,sub,0,dtin]

            # dtin=t3-tin-1
            # sign=1
            # if dtin<0:
                # dtin+=_map.MaxTauBin
                # sign*=-1
            # G1+=0.5*sign*G.Data[UP,sub,UP,sub,0,dtin]

            G1=sign*G.Data[UP,sub,UP,sub,0,dtin]

            for tout in range(_map.MaxTauBin):
                dtout=tout-t3-1
                sign=1
                if dtout<0:
                    dtout+=_map.MaxTauBin
                    sign*=-1
                # G2=0.5*sign*G.Data[UP,sub,UP,sub,0,dtout]

                # dtout=tout-t3
                # sign=1
                # if dtout<0:
                    # dtout+=_map.MaxTauBin
                    # sign*=-1
                # G2+=0.5*sign*G.Data[UP,sub,UP,sub,0,dtout]

                G2=sign*G.Data[UP,sub,UP,sub,0,dtout]

                GG=G1*G2
                
                for r1,r2,V in Neighbors:
                    GammaGNew[UP,r1,tout,tin]+=GG*(V[spinUP,spinUP]*GammaG[UP,r2,t3,t3]+V[spinUP,spinDOWN]*GammaG[DOWN,r2,t3,t3])
                    GammaGNew[DOWN,r1,tout,tin]+=GG*(V[spinDOWN,spinUP]*GammaG[UP,r2,t3,t3]+V[spinDOWN,spinDOWN]*GammaG[DOWN,r2,t3,t3])

    GammaGNew*=_map.Beta/_map.MaxTauBin
    return GammaGNew
