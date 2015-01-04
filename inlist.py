'''This is the input file of all jobs. 
   You have to add new job objects to TO_DO list
   if you want to run simulation.'''
import job_class as job
CPU = 4
SLEEP = 5    #check job status for every SLEEP seconds
TO_DO = []

#common dictionary for all jobs
beta=0.5
com_dict={
    "L" :   [8,8],
    "InitialBeta" :  beta,
    "DeltaBeta" :  0.00,
    "FinalBeta" :  beta,
    "Order" :  1,
    "Model" : "J1J2",
    #"Lattice" : "Square",
    "Lattice" : "Checkboard",
    #"Square", "Checkboard", "Honeycomb", "Square"
    "NSublat" : 2,
    "MaxTauBin" : 64,
    "Interaction" : [1.0,0.0],
    "ExternalField": 0.0,
}

# monte carlo job defintion
mc_dict={
    "__Execute" : "./gamma3.exe",
    "__Duplicate" : 1,
    "__IsCluster" : False,
    "__AutoRun" : True,
    "DoesLoad" : False,
    #Start from order 0, so that OrderReWeight has Order+1 elements
    "OrderReWeight" : [1.0, 1.0],
    "Sample" :  500000,
    "Sweep" : 10,
    "Toss" : 10000,
    "WormSpaceReweight" : 0.500
    }
mc_dict.update(com_dict)
TO_DO.append(job.JobMonteCarlo(mc_dict))

# self consist loop job definition
sc_dict={
    "__Execute" : ["python", "./run_loop.py"],
    "__Duplicate" : 1,
    "__IsCluster" : False,
    "__AutoRun" : False, 
    "DoesLoad" : True,
    "StartFromBare" : True,
    "OrderAccepted": 1,
    "ErrorThreshold": 0.5,
    "SleepTime": 300
    }
sc_dict.update(com_dict)
TO_DO.append(job.JobConsistLoop(sc_dict))

#diagram counter job definition
#diagcount_dict{
    #"__Execute" : "./gamma3.exe",
    #"__Duplicate" : 1,
    #"__IsCluster" : False,
    #"__AutoRun" : True,
    #"DoesLoad" : False,
    #Start from order 0, so that OrderReWeight has Order+1 elements
    #"OrderReWeight" : [1.0, 1.0, 3.0,4.0,1.0],
    #"Sample" :  5000000,
    #"Sweep" : 10,
    #"Toss" : 10000,
    #"WormSpaceReweight" : 0.100
#}

if __name__ == "__main__":
    for e in TO_DO:
        print e
        print e.to_dict(1)+"\n"

