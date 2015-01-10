'''This is the input file of all jobs. 
   You have to add new job objects to TO_DO list
   if you want to run simulation.'''
CPU = 4
SLEEP = 1    #check job status for every SLEEP seconds
#common dictionary for all jobs
beta=0.1
Common={
"Tau": {
    "MaxTauBin" : 32,
    "Beta": beta,
    "DeltaBeta" :  0.00,
    "FinalBeta" :  beta,
    },
"Lattice":  {
    "Name": "Square",
    "NSublat": 1,
    "L": [8,8]
    },
"Model": {
    "Name": "J1J2",
    "Interaction": [1.0,0.0],
    "ExternalField": [0.0]
    #ExternalField on Sublattice A and B
    }
}
# monte carlo job defintion
MonteCarlo={
"Control": {
    "__Execute" : "./simulator.exe",
    "__Duplicate" : 1,
    "__IsCluster" : False,
    "__AutoRun" : False,
    "__KeepCPUBusy": True,
    },
"Job": {"DoesLoad" : False},
"Markov": {
    "Order": 3,
    #Start from order 0, so that OrderReWeight has Order+1 elements
    "OrderReWeight" : [1.0, 1.0, 50.0, 500.0],
    "Sample" : 1000000,
    "Sweep" : 50,
    "Toss" : 10000,
    "WormSpaceReweight" : 0.05
    }
}
Dyson={
"Control": {
    "__Execute" : ["python", "./calculator/main.py"],
    "__Duplicate" : 1,
    "__IsCluster" : False,
    "__AutoRun" : True, 
    "__KeepCPUBusy": False,
    },
"Job": {"StartFromBare" : True},
"Dyson": {
    "Order": 1,
    "OrderAccepted": 1,
    "ErrorThreshold": 0.5,
    "SleepTime": 300
    }
}
import job_class as job
TO_DO = []
MonteCarlo.update(Common)
TO_DO.append(job.JobMonteCarlo(MonteCarlo))
Dyson.update(Common)
TO_DO.append(job.JobDyson(Dyson))
