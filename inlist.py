'''This is the input file of all jobs. 
   You have to add new job objects to TO_DO list
   if you want to run simulation.'''
import job_class as job
CPU = 4
SLEEP = 5    #check job status for every SLEEP seconds
TO_DO = []

#common dictionary for all jobs
com_dict={
    "L" :   [8,8],
    "Jcp" :  1.0,
    "iniBeta" :  1.00,
    "dBeta" :  0.00,
    "finalBeta" :  1.00,
    "Order" :  4,
    }

readfile="{0:4.2f}_{1}_coll".format(com_dict["finalBeta"],com_dict["Order"])
print readfile

# monte carlo job defintion
mc_dict={
    "__Execute" : "./gamma3.exe",
    "__Duplicate" : 4,
    "__IsCluster" : False,
    "__AutoRun" : True,
    "DoesLoad" : False,
    "Reweight" : [1.5, 1.0, 3.0,4.0],
    "ReadFile" : readfile,
    "Sample" :  5000000,
    "Sweep" : 10,
    "Toss" : 10000,
    "WormSpaceReweight" : 0.100
    }
mc_dict.update(com_dict)
TO_DO.append(job.JobMonteCarlo(mc_dict))

# self consist loop job definition
sc_dict={
    "__Execute" : ["python", "./run_loop.py"],
    "__Duplicate" : 0,
    "__IsCluster" : False,
    "__AutoRun" : False, 
    "DoesLoad" : True,
    "ReadFile" : readfile,
    }
sc_dict.update(com_dict)
TO_DO.append(job.JobConsistLoop(sc_dict))

# self consist loop job to initialize the simulation
sc_ini_dict={
    "__Execute" : ["python", "./run_loop.py"],
    "__Duplicate" : 0,
    "__IsCluster" : False,
    "__AutoRun" : False, 
    "DoesLoad" : False,
    "ReadFile" : readfile,
    }
sc_ini_dict.update(com_dict)

TO_DO.append(job.JobConsistLoop(sc_ini_dict))

# output loop job definition
ol_dict={
    "__Execute" : ["python", "./run_loop.py"],
    "__Duplicate" : 0,
    "__IsCluster" : False,
    #"__AutoRun" : True,
    "__AutoRun" : False,
    "DoesLoad" : True,
    "ReadFile" : readfile,
    }
ol_dict.update(com_dict)
TO_DO.append(job.JobOutputLoop(ol_dict))

# output numerical integration job definition
ni_dict={
    "__Execute" : ["./gamma3.exe"],
    "__Duplicate" : 0,
    "__IsCluster" : False,
    "__AutoRun" : False,
    "DoesLoad" : False,
    }
ni_dict.update(com_dict)
TO_DO.append(job.JobIntegration(ni_dict))


# output order job definition
oo_dict={
    "__Execute" : ["./gamma3.exe"],
    "__Duplicate" : 0,
    "__IsCluster" : False,
    "__AutoRun" : False,
    "DoesLoad" : True,
    "ReadFile" : readfile,
    }
oo_dict.update(com_dict)
TO_DO.append(job.JobOutputOrder(oo_dict))

# debug job definition
bg_dict={
    "__Execute" : ["./gamma3.exe"],
    "__Duplicate" : 0,
    "__IsCluster" : False,
    "__AutoRun" : False,
    "DoesLoad" : False,
    }
bg_dict.update(com_dict)
TO_DO.append(job.JobDebug(bg_dict))

if __name__ == "__main__":
    for e in TO_DO:
        print e
        print e.to_string(1)+"\n"

