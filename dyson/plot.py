#!/usr/bin/env python
from logger import *
import lattice as lat
import numpy as np
import weight, os, matplotlib
if "DISPLAY" not in os.environ:
    log.info("no DISPLAY detected, switch to Agg backend!")
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def PlotArray(array, Beta, Name, DoesSave=True):
    x=np.linspace(0, Beta, len(array))
    plt.figure()
    plt.plot(x,array,'-')
    if DoesSave:
        plt.savefig("{0}.jpg".format(Name))
    else:
        plt.show()
    plt.close()

def PlotTime(Name, weight, SpinIn, SubIn, SpinOut, SubOut, Vol, DoesSave=True):
    x=np.linspace(0, weight.Map.Beta, weight.Map.MaxTauBin)
    plt.figure(1)
    plt.suptitle("{0}, Spin ({1},{2}), Sublat ({3},{4}), Coordinates: {5}".format(Name,
        SpinIn, SpinOut, SubIn, SubOut, weight.Map.IndexToCoordi(Vol)))
    ax1=plt.subplot(121)
    ax1.plot(x,weight.Data[SpinIn, SubIn, SpinOut, SubOut, Vol,:].real,'-')
    ax1.set_xlabel("$\\tau$")
    ax1.set_ylabel("{0}.real".format(Name))
    ax2=plt.subplot(122)
    ax2.plot(x,weight.Data[SpinIn, SubIn, SpinOut, SubOut, Vol,:].imag,'-')
    ax2.set_xlabel("$\\tau$")
    ax2.set_ylabel("{0}.imag".format(Name))
    if DoesSave:
        plt.savefig("{0}.pdf".format(Name))
    else:
        plt.show()
    plt.close()

def PlotSpatial(weight, lattice, SpinIn, SpinOut, DoesSave=True):
    Omega=0
    OriginalSubLat=0
    coord3=0
    weight.FFT("R","T")
    x=[]
    y=[]
    z=[]
    points, _=lattice.GetSitesList()
    for vec, coord, sub in points:
        if lattice.Dim==2 or (lattice.Dim==3 and coord[2]==0):
            x.append(vec[0])
            y.append(vec[1])
            if coord[0]==coord[1]==0 and sub==OriginalSubLat:
                z.append(0.0)
            else:
                z.append(weight.Data[SpinIn,OriginalSubLat,SpinOut,sub,
                    weight.Map.CoordiIndex(coord),Omega].real)
    log.info("Max:{0}, Min: {1}".format(max(z), min(z)))
    plt.figure()
    plt.scatter(x,y,c=z, s=10, edgecolor="black", linewidth=0)
    c = plt.colorbar(orientation='horizontal', shrink=0.8, ticks=np.linspace(min(z),max(z),4))
    c.set_label("magnitude")
    plt.axis('equal')
    if DoesSave:
        plt.savefig("spatial_sub{0}.pdf".format(OriginalSubLat))
    else:
        plt.show()
    plt.close()

#def Exp(t, a, b):
    #return a*np.exp(-b*np.array(t))
def PlotWeightvsR(Name, weight, lattice, SpinIn, SpinOut, Tau=0, DoesSave=True):
    Omega=0
    OriginalSubLat=0
    coord3=0
    #weight.FFT("R","W")
    weight.FFT("R","T")
    x=[]
    y=[]
    z=[]
    points, _=lattice.GetSitesList(HasOffset=True)
    Norm = weight.Data[SpinIn, OriginalSubLat, SpinOut, OriginalSubLat, 0, Omega]
    for vec, coord, sub in points:
        if (all(v == 0 for v in coord) and sub==OriginalSubLat):
            Origin=np.array(vec)
    for vec, coord, sub in points:
        if not (all(v == 0 for v in coord) and sub==OriginalSubLat) and all(v<l/2 for v,l in zip(coord, lattice.L)):
        #if not (all(v == 0 for v in coord) and sub==OriginalSubLat):
            x.append(np.linalg.norm(np.array(vec)-Origin))
            #y.append(abs(weight.Data[SpinIn, OriginalSubLat, SpinOut, sub,
                #weight.Map.CoordiIndex(coord), Omega])/Norm)
            y.append(abs(weight.Data[SpinIn, OriginalSubLat, SpinOut, sub,
                weight.Map.CoordiIndex(coord), Omega]/Norm))
    #sort x,y according to the distance in x
    x,y = (list(x) for x in zip(*sorted(zip(x, y), key=lambda pair: pair[0])))
    #fitting
    #fitParams, fitCovariances = curve_fit(Exp, x, y)
    plt.figure()
    plt.plot(x,y, "o")
    #plt.plot(x, Exp(x, fitParams[0], fitParams[1]), '-',
            #label="fit with ${0}exp(-R/{1}a)$".format(fitParams[0], 1.0/fitParams[1]))
    plt.yscale("log")
    plt.xlabel("$R/a$")
    plt.ylabel("$|{0}(\omega={1})|$".format(Name,Omega))
    #plt.legend()
    if DoesSave:
        plt.savefig("WeightvsR.pdf")
    else:
        plt.show()
    plt.close()

def PlotChiAlongPath(Chi, lat, DoesSave=True):
    omega=0
    map=Chi.Map
    Chi.FFT("R", "W")
    #fig=plt.figure(figsize=(20, 10))
    fig=plt.figure()
    for BZcenter in lat.IndependtBZCenter:
        x=[]
        KList=[]
        offset=0
        ticks=[0]
        for i in range(0, len(lat.Path)-1):
            start, end=np.array(lat.Path[i]),np.array(lat.Path[i+1])
            for k in lat.GetKVecAlongPath(start, end, BZcenter):
                pos=offset+np.linalg.norm(k-np.array(BZcenter)-start)
                x.append(pos)
                KList.append(k)
            offset+=np.linalg.norm(end-start)
            ticks.append(offset)
        _, y=lat.FourierTransformation(Chi.Data[0,:,0,:,:,omega]*map.Beta/map.MaxTauBin, KList, "Real")
        y=[e.real for e in y]
        BZstr=["{:.3f}".format(e) for e in BZcenter]
        #x obtained previously may from big to small, so we have to reorder x here
        x,y=zip(*sorted(zip(x, y)))
        plt.plot(x,y,'o-', label="BZ:({0})".format(",".join(BZstr)))
    plt.legend(loc='best', fancybox=True, framealpha=0.5)
    plt.xticks(ticks, lat.PathName)
    if DoesSave:
        plt.savefig("chi_1D.pdf")
    else:
        plt.show()
    plt.close()

def PlotChi_2D(Chi, lat, DoesSave=True):
    omega=0
    map=Chi.Map
    Chi.FFT("R", "W")

    if lat.Name=="Pyrochlore":
        #####Pyrochlore
        KList_hhl=[]
        KList_hl0=[]
        for i in range(-lat.L[0]*4, lat.L[0]*4+1):
            for j in range(-lat.L[1]*4, lat.L[1]*4+1):
                for k in range(-lat.L[2]*4, lat.L[2]*4+1):
                    kpoint = i*lat.ReciprocalLatVec[0]+j*lat.ReciprocalLatVec[1]+ \
                            k*lat.ReciprocalLatVec[2]
                    if np.abs(kpoint[0]-kpoint[1])<1e-5:
                        KList_hhl.append((i,j,k))
                    if np.abs(kpoint[2])<1e-5:
                        KList_hl0.append((i,j,k))

        bound=[[-40,40],[-40,40]]
        ######hhl
        k_hhl, ChiK_hhl=lat.FourierTransformation(Chi.Data[0,:,0,:,:,omega]*map.Beta/map.MaxTauBin, \
                KList_hhl, "Integer", bound=bound)
        ChiK_hhl=[e.real for e in ChiK_hhl]

        x_hhl=[]
        y_hhl=[]
        for e in k_hhl:
            x_hhl.append(np.sqrt(2.0)*e[0])
            y_hhl.append(e[2])

        ######hl0
        k_hl0, ChiK_hl0=lat.FourierTransformation(Chi.Data[0,:,0,:,:,omega]*map.Beta/map.MaxTauBin, \
                KList_hl0, "Integer", bound=bound)
        ChiK_hl0=[e.real for e in ChiK_hl0]
        x_hl0=[]
        y_hl0=[]

        for e in k_hl0:
            x_hl0.append(e[0])
            y_hl0.append(e[1])

        plt.figure(1)
        ax1=plt.subplot(121,aspect='equal')
        plt.scatter(x_hhl,y_hhl,c=ChiK_hhl, s=29, edgecolor="black", linewidth=0)
        plt.xlabel("Direction [hh0]")
        plt.ylabel("Direction [00l]")
        plt.xlim(-30, 30)
        plt.ylim(-30, 30)
        label=np.linspace(min(ChiK_hhl),max(ChiK_hhl), 4)

        PI2=2*np.pi
        sqrt2 = np.sqrt(2.0)
        xlist = PI2*sqrt2*np.array([-0.75,-0.25, 0.25, 0.75, 0.25,-0.25,-0.75])
        ylist = PI2*np.array([          0,    1,    1,    0,   -1,   -1,    0])
        lc="white"
        plt.plot(xlist, ylist, color=lc)
        plt.plot(xlist, ylist+2*PI2, color=lc)
        plt.plot(xlist, ylist-2*PI2, color=lc)
        plt.plot(xlist+sqrt2*PI2, ylist+1*PI2, color=lc)
        plt.plot(xlist-sqrt2*PI2, ylist+1*PI2, color=lc)
        plt.plot(xlist+sqrt2*PI2, ylist-1*PI2, color=lc)
        plt.plot(xlist-sqrt2*PI2, ylist-1*PI2, color=lc)
        c = plt.colorbar(orientation='horizontal', shrink=0.8, ticks=label)
        c.set_label("magnitude")

        ax2=plt.subplot(122,aspect='equal')
        plt.scatter(x_hl0,y_hl0,c=ChiK_hl0, s=18, edgecolor="black", linewidth=0)
        plt.xlabel("Direction [h00]")
        plt.ylabel("Direction [0l0]")
        label=np.linspace(min(ChiK_hl0),max(ChiK_hl0), 4)
        plt.xlim(-40, 40)
        plt.ylim(-40, 40)
        plt.tight_layout()

        xlist = PI2*np.array([-1.0,-0.5, 0.5, 1.0, 1.0, 0.5,-0.5,-1.0,-1.0])
        ylist = PI2*np.array([ 0.5, 1.0, 1.0, 0.5,-0.5,-1.0,-1.0,-0.5, 0.5])
        plt.plot(xlist, ylist, color=lc)
        plt.plot(xlist+2*PI2, ylist, color=lc)
        plt.plot(xlist, ylist+2*PI2, color=lc)
        plt.plot(xlist+2*PI2, ylist+2*PI2, color=lc)
        c = plt.colorbar(orientation='horizontal', shrink=0.8, ticks=label)
        c.set_label("magnitude")

    elif lat.Name in ["3DCheckerboard", "Cubic"]:
        ####3D Checkerboard
        KList_hl0=[]
        KList_hhl=[]

        for i in range(-2*lat.L[0]+1, 2*lat.L[0]):
            for j in range(-2*lat.L[1]+1, 2*lat.L[1]):
                KList_hl0.append((i,j,0))

        k_hl0, ChiK_hl0=lat.FourierTransformation(Chi.Data[0,:,0,:,:,omega]*map.Beta/map.MaxTauBin,
                KList_hl0, "Integer")
        ChiK_hl0=[e.real for e in ChiK_hl0]
        x_hl0=[]
        y_hl0=[]
        for e in k_hl0:
            x_hl0.append(e[0])
            y_hl0.append(e[1])

        plt.figure(1)
        plt.scatter(x_hl0,y_hl0,c=ChiK_hl0, s=10, edgecolor="black", linewidth=0)
        c = plt.colorbar(orientation='horizontal')
        c.set_label("magnitude")
        plt.axis('equal')

    elif lat.Dim==2:
        KList=[]
        for i in range(-2*lat.L[0], 2*lat.L[0]+1):
            for j in range(-2*lat.L[1], 2*lat.L[1]+1):
                KList.append((i,j))
        k, ChiK=lat.FourierTransformation(Chi.Data[0,:,0,:,:,omega]*map.Beta/map.MaxTauBin,
                KList, "Integer", bound=[[-20,20], [-20,20]])
        ChiK=[e.real for e in ChiK]
        k=np.array(k)
        plt.figure()
        plt.scatter(k[:, 0],k[:, 1],c=ChiK, s=6, edgecolor="black", linewidth=0)
        c = plt.colorbar(orientation='horizontal')
        c.set_label("magnitude")
        plt.axis('equal')
        Ktemp=[(-3,0),(3,0),(0,3),(0,-3)]
        k, ChiK=lat.FourierTransformation(Chi.Data[0,:,0,:,:,omega]*map.Beta/map.MaxTauBin,
                Ktemp, "Integer")
    else:
        log.warn("Lattice PlotChi_2D not implemented yet!")

    if DoesSave:
        plt.savefig("chiK_{0}.pdf".format(lat.Name))
    else:
        plt.show()
    plt.close()
    log.info("Plotting done!")


if __name__=="__main__":
    import weight
    import IO

    WeightPara={"NSublat": 4, "L":[8,8,8],
            "Beta": 6.0, "MaxTauBin":64}
    Map=weight.IndexMap(**WeightPara)
    l=lat.Lattice("Pyrochlore", Map)

    Dict = IO.LoadBigDict("Weight")["Chi"]
    Chi=weight.Weight("SmoothT", Map, "NoSpin", "Symmetric","R","T").FromDict(Dict)

    PlotChiAlongPath(Chi, l)
    PlotChi_2D(Chi, l)
    PlotWeightvsR("\chi", Chi,l,0,0)

