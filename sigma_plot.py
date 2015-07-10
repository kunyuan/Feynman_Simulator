#!/usr/bin/env python
import read_data
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

Beta = 1.25
N =128

tau = np.arange(Beta/(2*N), Beta+Beta/(2*N), Beta/N)
target=["Sigma"]

BoldSigma=[]
BoldSigma.append(read_data.read_array("../L16_1.25_2/1.25_quantities.dat",target))
BoldSigma.append(read_data.read_array("../data_Nikolay/Sigma_Order2_Nikolay.dat",target))

fig = plt.figure()
ax = plt.subplot(111)

for key in target:
    ax.plot(tau, BoldSigma[0][key][0].real/16.0, label="{0} Yuan".format(key))
    #ax.plot(tau, BoldSigma[0][key][0].imag/16.0, label="{0} Yuan".format(key))
    ax.plot(tau, -BoldSigma[1][key][0].real, label="{0} Nikolay".format(key))
    #ax.plot(tau, BoldSigma[1][key][0].imag, label="{0} Nikolay".format(key))

ax.legend()

plt.xlabel("tau")
plt.ylabel("Sigma")

plt.savefig("Beta1.25_L16_Sigma_full.pdf")
plt.show()
