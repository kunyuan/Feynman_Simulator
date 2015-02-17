'''This is the input file of all jobs. 
   You have to add new job objects to TO_DO list
   if you want to run simulation.'''
CPU = 4
SLEEP = 1    #check job status for every SLEEP seconds
#common dictionary for all jobs
# monte carlo job defintion
MonteCarlo={
"Control": {
    "__Execute" : "./simulator.exe",
    "__Duplicate" :  3,
    "__IsCluster" : False, 
    "__AutoRun" : True,
    "__KeepCPUBusy": True,
    },
"Job": {
    "DoesLoad" : False,
    "Sample" : 100000000  ##0.8 min for 1000000(*1000) Samples in MC
    }
}

Dyson={
"Control": {
    "__Execute" : ["python", "./dyson/main.py"],
    "__Duplicate" : 1,
    "__IsCluster" : MonteCarlo["Control"]["__IsCluster"],
    "__AutoRun" : True, 
    "__KeepCPUBusy": False,
    },
"Job": {
    "DysonOnly": MonteCarlo["Control"]["__Duplicate"]==0
    #"DysonOnly": False
    }
}

<<<<<<< HEAD
Beta=4.0
Order=4
=======
Beta=2.0
Order=1
>>>>>>> 877f720202e931f7a6d27e5753d6ba8d8f15d647
Common={
"Tau": {
    "MaxTauBin" : 64,
    "Beta": Beta,
    },
"Lattice":  {
    #"Name": "Square",
    #"NSublat": 1,
    #"L": [16,16],
    #"Name": "Honeycomb",
    #"NSublat": 2,
    #"L": [16,16],
    #"Name": "Kagome",
    #"NSublat": 3,
    #"L": [16,16],
    #"Name": "Cubic",
    #"NSublat": 1,
    #"L": [8,8,8],
    "Name": "Pyrochlore",
    "NSublat": 4,
    "L": [8,8,8]
    #"Name": "Checkboard",
    #"NSublat": 2,
    #"L": [8,8]
    #"Name": "3DCheckerboard",
    #"NSublat": 2,
    #"L": [16,16,16]
    #"L": [8,8,8]
    },
"Model": {
    "Name": "J1J2",
    #"Description": ["ImW",],
    "Interaction": [1.0,0.0],
    "ExternalField": [0.0, 0.0, 0.0, 0.0]
    #ExternalField on Sublattice A and B
    },
}

MonteCarlo["Markov"]={
    "Order": Order,
    #Start from order 0, so that OrderReWeight has Order+1 elements
    "Sweep" : 10,
    "Toss" : 1000,
    "OrderReWeight" : [100.0, 0.5, 1.0, 0.1, 0.05, 0.05, 0.01, 0.005],
    "WormSpaceReweight" : 0.05,
    "PolarReweight" : 2.0,
    "OrderTimeRatio" : [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    #"Timer": {
        #"PrinterTimer": 300,
        #"DiskWriterTimer": 300,
        #"MessageTimer": 310,
        #"ReweightTimer": 600
        #},
    "Timer": {
        "PrinterTimer": 60,
        "DiskWriterTimer": 60,
        "MessageTimer": 60,
        "ReweightTimer":60
        },
    }

Dyson["Dyson"]={
    "OrderAccepted": {"Sigma":1, "Polar":1},
    "ErrorThreshold": 0.1,
    #"SleepTime": 300,
    "SleepTime": 40,
    "Annealing": {
        "DeltaField": [0.0, 0.0, 0.0, 0.0],
        "Interval": [-0.1, -0.1, -0.1, -0.1]
        }
    }

import job_class as job
TO_DO = []
MonteCarlo.update(Common)
TO_DO.append(job.JobMonteCarlo(MonteCarlo))
Dyson.update(Common)
TO_DO.append(job.JobDyson(Dyson))
