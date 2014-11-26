'''This is the input file of all jobs. 
   You have to add new job objects to TO_DO list
   if you want to run simulation.'''
import job_class as job
import model
CPU = 4
SLEEP = 5    #check job status for every SLEEP seconds
TO_DO = []

#common dictionary for all jobs
beta=1.0
com_dict={
    "L" :   [8,8],
    "initialBeta" :  beta,
    "deltaBeta" :  0.00,
    "finalBeta" :  beta,
    "Order" :  4,
    }

hamiltonian_dict={
    "Interaction" : [1.0,0.5],
    "ExternalField": 0.0,
}
hamiltonian_dict.update(model.J1J2(beta))

# monte carlo job defintion
mc_dict={
    "__Execute" : "./gamma3.exe",
    "__Duplicate" : 1,
    "__IsCluster" : False,
    "__AutoRun" : True,
    "DoesLoad" : False,
    "StartFromBare" : True,
    
    "OrderReweight" : [1.5, 1.0, 3.0,4.0],
    "Sample" :  5000000,
    "Sweep" : 10,
    "Toss" : 10000,
    "WormSpaceReweight" : 0.100
    }
mc_dict.update(com_dict)
mc_dict.update(hamiltonian_dict)
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
sc_dict.update(hamiltonian_dict)
TO_DO.append(job.JobConsistLoop(sc_dict))

if __name__ == "__main__":
    for e in TO_DO:
        print e
        print e.to_string(1)+"\n"

