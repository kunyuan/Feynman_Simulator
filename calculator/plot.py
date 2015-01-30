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
    omega=1
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
        plt.show()
    plt.clf()
    plt.cla()
    plt.close()

    #KList=[]
    #for i in range(-2*lat.L[0], 2*lat.L[0]):
        #for j in range(-2*lat.L[1], 2*lat.L[1]):
            #KList.append((i,j,0))
    #k, ChiK=lat.FourierTransformation_RealSpace(Chi.Data[0,0,0,:,:,omega]*map.Beta/map.MaxTauBin, KList, "Integer")
    #x=[]
    #y=[]
    #for e in k:
        #x.append(e[0])
        #y.append(e[1])
    #plt.scatter(x,y,c=ChiK)
    #c = plt.colorbar(orientation='horizontal')
    #c.set_label("magnitude")
    #plt.axis('equal')
    #if DoesSave:
        #plt.savefig("chiK.jpg")
    #else:
        #plt.show()
    #plt.clf()
    #plt.cla()
    #plt.close()
    Chi.FFT("R", "T")


