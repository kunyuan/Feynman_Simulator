import sys
import IO
import matplotlib.pyplot as plt

dataMC=IO.LoadBigDict("data/2_statistics")['Histogram']['Sigma']['Histogram']
sigmaMC = dataMC['WeightAccu'][0][0][0][0]*dataMC['Norm']/dataMC['NormAccu']
plt.plot(sigmaMC[1:])

dataDyson=IO.LoadBigDict("data/Weight")['Sigma']
plt.plot(dataDyson['SmoothT'][0][0][0])
plt.show()
