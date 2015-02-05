import lattice as lat
import numpy as np
import weight
#from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from logger import *

def PlotArray(array, Beta, Name, DoesSave=True):
    x=np.linspace(0, Beta, len(array))
    plt.plot(x,array,'-')
    if DoesSave:
        plt.savefig("{0}.jpg".format(Name))
    else:
        plt.show()
    plt.clf()
    plt.cla()
    plt.close()

def PlotTime(weight, SpinIn, SubIn, SpinOut, SubOut, Vol, DoesSave=True):
    x=np.linspace(0, weight.Map.Beta, weight.Map.MaxTauBin)
    plt.plot(x,weight.Data[SpinIn, SubIn, SpinOut, SubOut, Vol,:],'-')
    if DoesSave:
        plt.savefig("time.jpg")
    else:
        plt.show()
    plt.clf()
    plt.cla()
    plt.close()

def PlotSpatial(weight, lattice, SpinIn, SpinOut, Tau=0, DoesSave=True):
    if lattice.Dim>=2:
        log.warning("Can not plot for Dimension>3!")
        return
    map=weight.Map
    x=[]
    y=[]
    z=[]
    points,_=lattice.GetSitesList()
    for vec, coord, sub in points:
        x.append(vec[0])
        y.append(vec[1])
        if coord[0]==coord[1]==sub==0:
            z.append(0.0)
        else:
            z.append(weight.Data[SpinIn,0,SpinOut,sub,map.CoordiIndex(coord),map.TauIndex(0.0,Tau)].real)
        #print z.real
    plt.scatter(x,y,s=100,c=z)
    c = plt.colorbar(orientation='horizontal')
    c.set_label("magnitude")
    plt.axis('equal')
    if DoesSave:
        plt.savefig("spatial.jpg")
    else:
        plt.show()
    plt.clf()
    plt.cla()
    plt.close()

def PlotChi(Chi, lat, DoesSave=True):
    omega=0
    map=Chi.Map
    x=[]
    KList=[]
    offset=0
    ticks=[0]
    Chi.FFT("R", "W")
    for i in range(0, len(lat.Path)-1):
        start=np.array(lat.Path[i])
        end=np.array(lat.Path[i+1])
        N=lat.PathNum[i]+1
        path=np.zeros((map.Dim, N))
        dpath=np.zeros((map.Dim, N))
        for j in range(0, map.Dim):
            p=np.linspace(start[j],end[j],N)
            path[j,:]=p
            dpath[j,:]=p-start[j]
        for k in range(0,N):
            x.append(offset+np.linalg.norm(dpath[:,k]))
            KList.append(path[:,k])
        offset+=np.linalg.norm(end-start)
        ticks.append(offset)
    k, y=lat.FourierTransformation_RealSpace(Chi.Data[0,0,0,:,:,omega]*map.Beta/map.MaxTauBin, KList, "Real")
    plt.plot(x,y,'o')
    plt.xticks(ticks, lat.PathName)
    if DoesSave:
        plt.savefig("chi.jpg")
    else:
        print "show plot"
        plt.show()
    plt.clf()
    plt.cla()
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
                #KList_hl0.append((i*2*np.pi/lat.L[0],j*2*np.pi/lat.L[1],0))

        ######hhl
        k_hhl, ChiK_hhl=lat.FourierTransformation_RealSpace(Chi.Data[0,0,0,:,:,omega]*map.Beta/map.MaxTauBin, KList_hhl, "Integer")

        x_hhl=[]
        y_hhl=[]
        for e in k_hhl:
            kpoint  = e[0] * lat.ReciprocalLatVec[0]/2/np.pi
            kpoint += e[1] * lat.ReciprocalLatVec[1]/2/np.pi
            kpoint += e[2] * lat.ReciprocalLatVec[2]/2/np.pi
            x_hhl.append(np.sqrt(2.0)*kpoint[0])
            y_hhl.append(kpoint[2])

        ######hl0
        k_hl0, ChiK_hl0=lat.FourierTransformation_RealSpace(Chi.Data[0,0,0,:,:,omega]*map.Beta/map.MaxTauBin, KList_hl0, "Integer")
        #k_hl0, ChiK_hl0=lat.FourierTransformation_RealSpace(Chi.Data[0,0,0,:,:,omega]*map.Beta/map.MaxTauBin, KList_hl0, "Real")
        x_hl0=[]
        y_hl0=[]

        for e in k_hl0:
            kpoint  = e[0] * lat.ReciprocalLatVec[0]/2/np.pi
            kpoint += e[1] * lat.ReciprocalLatVec[1]/2/np.pi
            kpoint += e[2] * lat.ReciprocalLatVec[2]/2/np.pi
            x_hl0.append(kpoint[0])
            y_hl0.append(kpoint[1])

            #x_hl0.append(e[0])
            #y_hl0.append(e[1])

        plt.figure(1)
        plt.subplot(121)
        plt.scatter(x_hhl,y_hhl,c=ChiK_hhl, s=10, edgecolor="black", linewidth=0)
        c = plt.colorbar(orientation='horizontal')
        c.set_label("magnitude")
        #plt.axis('equal')

        plt.subplot(122)
        plt.scatter(x_hl0,y_hl0,c=ChiK_hl0, s=10, edgecolor="black", linewidth=0)
        c = plt.colorbar(orientation='horizontal')
        c.set_label("magnitude")
        plt.axis('equal')

        if DoesSave:
            plt.savefig("chiK_Pyrochlore.jpg")
        else:
            plt.show()
        plt.clf()
        plt.cla()
        plt.close()

    elif lat.Name=="Checkerboard":
        ####2D Checkerboard
        KList=[]
        for i in range(-2*lat.L[0], 2*lat.L[0]):
            for j in range(-2*lat.L[1], 2*lat.L[1]):
                KList.append((i,j,0))
        k, ChiK=lat.FourierTransformation_RealSpace(Chi.Data[0,0,0,:,:,omega]*map.Beta/map.MaxTauBin, KList, "Integer")
        x=[]
        y=[]
        for e in k:
            x.append(e[0])
            y.append(e[1])
        plt.scatter(x,y,c=ChiK)
        c = plt.colorbar(orientation='horizontal')
        c.set_label("magnitude")
        plt.axis('equal')
        if DoesSave:
            plt.savefig("chiK_2DChecker.jpg")
        else:
            plt.show()
        plt.clf()
        plt.cla()
        plt.close()

    elif lat.Name=="3DCheckerboard":
        ####3D Checkerboard
        KList_hl0=[]
        KList_hhl=[]

        for i in range(-2*lat.L[0], 2*lat.L[0]+1):
            for j in range(-2*lat.L[1], 2*lat.L[1]+1):
                #for k in range(-2*lat.L[1], 2*lat.L[1]+1):
                KList_hl0.append((i*2*np.pi/lat.L[0],0,j*2*np.pi/lat.L[2]))

        k, ChiK=lat.FourierTransformation_RealSpace(Chi.Data[0,0,0,:,:,omega]*map.Beta/map.MaxTauBin, KList_hl0, "Real")

        x=[]
        y=[]
        for e in k:
            x.append(e[0])
            y.append(e[2])

        plt.scatter(x,y,c=ChiK)
        c = plt.colorbar(orientation='horizontal')
        c.set_label("magnitude")
        plt.axis('equal')
        if DoesSave:
            plt.savefig("chiK_3DChecker.jpg")
        else:
            plt.show()
        plt.clf()
        plt.cla()
        plt.close()
    else:
        log.warn("Lattice PlotChi_2D not implemented yet!")
    log.info("Plotting done!")


