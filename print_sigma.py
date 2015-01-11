import sys
import IO
import read_data
import matplotlib.pyplot as plt

input = IO.LoadDict("infile/_in_MC_1")['Para']
MaxTauBin = input['Tau']['MaxTauBin']
Beta = input['Tau']['Beta']
Vol = input['Lattice']['L'][0]*input['Lattice']['L'][1]
SublatVol = input['Lattice']['NSublat']
normratio = Vol*(MaxTauBin/Beta)/Beta/Vol/SublatVol

dataMC=IO.LoadBigDict("1_MC_statis")['Sigma']['Histogram']['SmoothT']
sigmaMC={}
sigmaMC[0] = dataMC['WeightAccu'][0][0][0][0]*dataMC['Norm']/dataMC['NormAccu']
sigmaMC[1] =dataMC['WeightAccu'][1][0][0][0]*dataMC['Norm']/dataMC['NormAccu']
sigmaMC[2] =dataMC['WeightAccu'][2][0][0][0]*dataMC['Norm']/dataMC['NormAccu']

#plt.plot(normratio*sigmaMC[0].real)
plt.plot(normratio*sigmaMC[1].real)
plt.plot(normratio*sigmaMC[2].real)

target=["Sigma0", "Sigma1", "Sigma2", "Sigma"]
dataNikolay = read_data.read_array("Sigma_Order2_Nikolay.dat", target)
#plt.plot(-dataNikolay["Sigma0"][0].real)
plt.plot(-dataNikolay["Sigma1"][0].real)
plt.plot(-dataNikolay["Sigma2"][0].real)

#dataDyson=IO.LoadBigDict("Weight")['Polar']
#plt.plot(dataDyson['SmoothT'][0][0][0].real)

plt.show()
