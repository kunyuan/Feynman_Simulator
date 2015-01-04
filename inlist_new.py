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
"MaxTauBin" : 64,
"Beta": beta,
"L": [8,8],
"Model": {
    "Name": "J1J2",
    "Lattice":  {
        "Name": "Checkboard",
        "NSublat": 2
        }
    "Hamiltonian": {
        "Interaction": [1.0,0.0],
        "ExternalField": 0.0
        }
    }
}

# monte carlo job defintion
mc_dict={
    "__Execute" : "./gamma3.exe",
    "__Duplicate" : 1,
    "__IsCluster" : False,
    "__AutoRun" : True,
    "DoesLoad" : False,
    #Start from order 0, so that OrderReWeight has Order+1 elements
"MaxOrder": 1,
"OrderReWeight" : [1.0, 1.0],
"Sample" :  500000,
"Sweep" : 10,
"Toss" : 10000,
"WormSpaceReweight" : 0.500
}
mc_dict.update(com_dict)
TO_DO.append(job.JobMonteCarlo(mc_dict))

