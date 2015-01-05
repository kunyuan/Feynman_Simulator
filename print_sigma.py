import sys
import IO
import matplotlib.pyplot as plt

dataMC=IO.LoadBigDict("data/2_statistics")['Sigma']['Histogram']
sigmaMC =dataMC['WeightAccu'][0][0][0][0]*dataMC['Norm']/dataMC['NormAccu']
plt.plot(0.3*64*sigmaMC.real)

dataDyson=IO.LoadBigDict("data/Weight")['Sigma']
plt.plot(dataDyson['SmoothT'][0][0][0].real)
plt.show()
