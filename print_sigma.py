import numpy as np
import matplotlib.pyplot as plt
dataMC=np.load("data/2_statistics.npz")
sigmaMC = dataMC['Sigma.WeightAccu'][0][0][0][0]*dataMC['Sigma.Norm']/dataMC['Sigma.NormAccu']
plt.plot(sigmaMC[1:])
dataDyson=np.load("data/Weight.npz")
plt.plot(dataDyson['Sigma.SmoothT'][0][0][0])
plt.show()
