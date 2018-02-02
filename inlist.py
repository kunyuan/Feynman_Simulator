# monte carlo job defintion
MonteCarlo={
"Control": {
    "__Execute" : "./simulator.exe",
    "__Duplicate" :  0,
    # "__Duplicate" :  14,
    "__IsCluster" : False, 
    "__AutoRun" : True,
    },
"Job": {"Sample" : 100000000}  ##0.8 min for 1000000(*1000) Samples in MC
}
Dyson={
"Control": {
    "__Execute" : ["python", "./dyson/main.py"],
    "__Duplicate" : 1,
    "__IsCluster" : MonteCarlo["Control"]["__IsCluster"],
    # "__AutoRun" : MonteCarlo["Control"]["__AutoRun"], 
    "__AutoRun" : False,
    "__PBSCommand": "#PBS -l mem=5gb"
    },
"Job": {
    "DysonOnly": MonteCarlo["Control"]["__Duplicate"]==0,
    #"DysonOnly": False,
    "SumRule": False 
    # "SumRule":  True
    }
}

Beta=1.0
Order=3
MaxTauBin=64
L=8

Gamma3=False
# MaxTauBinTiny=MaxTauBin #the tau bin for large object like GammaW
MaxTauBinTiny=MaxTauBin/2 #the tau bin for large object like GammaW
BasisNum=16 #the number of basis

Common={
"Gamma3": Gamma3,
"Tau": {"MaxTauBin" : MaxTauBin, "Beta": Beta, "MaxTauBinTiny": MaxTauBinTiny, "BasisNum": BasisNum},
"Lattice":  {
    #1D lattice
    # "Name": "Chain", "NSublat":1,
    # "L": [2,]

    #2D lattice
    # "Name": "Square", "NSublat": 1,
    # "Name": "Checkerboard", "NSublat": 2,
    # "Name": "ValenceBond", "NSublat": 2,
    # "Name": "Honeycomb", "NSublat": 2,
    #"Name": "Kagome", "NSublat": 3,
    # "Name": "Triangular", "NSublat": 1,
    # "Name": "Assymetric_Triangular", "NSublat": 3,
    # "L": [L,L]

    #3D lattice
    #"Name": "Cubic", "NSublat": 1,
    #"Name": "3DCheckerboard", "NSublat": 2,
    "Name": "Pyrochlore", "NSublat": 4,
    "L": [L,L,L]
    },
"Model": {
    "Name": "J1J2",
    # "Name": "Kitaev",
    # "Name": "Heisenberg",
    #"Description": ["ImW",],
    "Interaction": [1.0, 0.0, 0.0, 0.0],
    "ExternalField": [ 0.0, 0.0, 0.0, 0.0],
    #ExternalField on Sublattice A and B
    "Symmetry": ["SU2", "ParticleHole"],
    }
}

MonteCarlo["Markov"]={
    "Order": Order, "Sweep" : 10, "Toss" : 1000,
    #Start from order 0, so that OrderReWeight has Order+1 elements
    "OrderReWeight" : [1.0, 0.3, 1.0, 0.1, 0.05, 0.05, 0.01, 0.005],
    "WormSpaceReweight" : 0.05,
    "PolarReweight" : 1.0,
    "runGamma3": Gamma3,
    "GammaGReweight" : 1.0,
    "GammaWReweight" : 1.0,
    "OrderTimeRatio" : [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    #"Timer": {
        #"PrinterTimer": 300,
        #"DiskWriterTimer": 300,
        #"MessageTimer": 310,
        #"ReweightTimer": 600
        #},
    "Timer": {
        "PrinterTimer": 30,
        "DiskWriterTimer": 30,
        "MessageTimer": 30,
        "ReweightTimer":1000
        },
    }

Dyson["Dyson"]={
    #"SleepTime": 0,
    "SleepTime": 60,
    "Gamma3": Gamma3,
    "OrderAccepted": {"Sigma":1, "Polar":1},
    "ErrorThreshold": 0.2,
    "Annealing": {
        #"DeltaField": [-0.5,  -0.5, 0.0, 0.0],
        #"Interval": [0.05, 0.05, -0.0, -0.0]
        "DeltaField": [-0.0,  -0.0, 0.0, 0.0],
        "Interval": [0.00, 0.00, -0.0, -0.0]
        }
    }

import job_class as job
'''This is the input file of all jobs. 
   You have to add new job objects to TO_DO list
   if you want to run simulation.'''
TO_DO = []
MonteCarlo.update(Common)
TO_DO.append(job.JobMonteCarlo(MonteCarlo))
Dyson.update(Common)
TO_DO.append(job.JobDyson(Dyson))
CPU = 16
SLEEP = 1    #check job status for every SLEEP seconds

import basis
svd=basis.SVDBasis(MaxTauBin,Beta, "Fermi")
svd.GenerateBasis(BasisNum)
svd.Save("FermiBasis.dat")
svd=basis.SVDBasis(MaxTauBin,Beta, "Bose")
svd.GenerateBasis(BasisNum)
svd.Save("BoseBasis.dat")
