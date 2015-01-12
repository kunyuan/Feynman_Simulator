import lattice as lat
import numpy as np
import weight
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from logger import *

def PlotTime(weight, SpinIn, SubIn, SpinOut, SubOut, Vol):
    x=np.linspace(0, weight.Map.Beta, weight.Map.MaxTauBin)
    plt.plot(x,weight.Data[SpinIn, SubIn, SpinOut, SubOut, Vol,:],'-')
    plt.savefig("time.jpg")
    plt.show()

def PlotSpatial(weight, lattice, SpinIn, SpinOut, Tau=0):
    Assert(lattice.Dim==2, "PlotSpatial only works for two dimensional system")
    map=weight.Map
    x=[]
    y=[]
    z=[]
    for vec, coord, sub in lattice.GetSitesList():
        x.append(vec[0])
        y.append(vec[1])
        if coord[0]==coord[1]==sub==0:
            z.append(0.0)
        else:
            z.append(weight.Data[SpinIn,0,SpinOut,sub,map.CoordiIndex(coord),map.TauIndex(0.0,Tau)].real)
        #print z.real
    print z
    plt.scatter(x,y,s=100,c=z)
    c = plt.colorbar(orientation='horizontal')
    c.set_label("magnitude")
    plt.savefig("spatial.jpg")
    plt.show()


    



