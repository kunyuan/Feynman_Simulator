import sys
import IO
import matplotlib.pyplot as plt

input = IO.LoadDict("data/infile/_in_MC_2")['Para']
MaxTauBin = input['Tau']['MaxTauBin']
Beta = input['Tau']['Beta']
Vol = input['Lattice']['L'][0]*input['Lattice']['L'][1]
SublatVol = input['Lattice']['NSublat']
normratio = (MaxTauBin/Beta)/Beta/Vol/SublatVol

dataMC=IO.LoadBigDict("data/2_statistics")['Sigma']['Histogram']
sigmaMC =dataMC['WeightAccu'][0][0][0][0]*dataMC['Norm']/dataMC['NormAccu']
plt.plot(normratio*sigmaMC.real)

dataDyson=IO.LoadBigDict("data/Weight")['Sigma']
plt.plot(dataDyson['SmoothT'][0][0][0].real)

plt.show()
