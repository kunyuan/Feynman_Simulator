''' This is a file to define the class of all jobs,
    You have to modify this file if you want to add new type of jobs.'''
import sys
import os
import random
#base class of all jobs
class Job:
    '''Class Job is the base class of all job objects.
       You have to add new subclass of Job if you want to
       create new type of job.'''
    def __init__(self, para):
        # interesting point here, python will actually call the 
        #__check_parameters__ of subclass,
        # So you don't have to call __check_parameters__ 
        #in subclass initialization again!
        if self.__check_parameters__(para) is False:
            print "Something is wrong with the inlist! Abandon!"
            sys.exit()
        self.duplicate = para.pop("__Duplicate")

        #take care of the list of paths
        execu = para.pop("__Execute")
        if type(execu) is str:
            execu = os.path.abspath(execu)
        else:
            for i in range(0, len(execu)):
                if os.path.isfile(execu[i]):
                    execu[i] = os.path.abspath(execu[i])
        #self.execute is the execute file str
        self.execute = execu
        self.is_cluster = para.pop("__IsCluster")
        self.auto_run = para.pop("__AutoRun")
        self.keep_cpu_busy = True
        self.pid = 0
        self.name = ""
        self.para = para

    def to_string(self, pid=0):
        '''output the corresponding string of the job class'''
        self.para["PID"] = pid
        self.para["WeightFile"]="Weight.npz"
        self.para["MessageFile"]="Message.txt"
        self.__set_model_specific__()
        return self.__formator__(self.para)

    def __formator__(self,para):
        import pprint
        return "Para="+pprint.pformat(para)

    def __check_parameters__(self, para):
        if para["__Execute"] is "":
            print "Please specify the executive file name!"
            return False

    def __set_model_specific__(self):
        PI=3.141592653589793238
        if self.para["Model"]=="J1J2":
            self.para["Hopping"]=[0.0,]
            mu=(0,PI/2.0/self.para["InitialBeta"])
            self.para["ChemicalPotential"]=[mu,mu]

class JobMonteCarlo(Job):
    '''job subclass for monte carlo jobs'''
    def __init__(self, para):
        Job.__init__(self, para)
        self.keep_cpu_busy = True
        self.para["Type"] = "MC"
        self.name = "MC"

    def __check_parameters__(self, para):
        if Job.__check_parameters__(self, para) is False:
            return False
        if type(para["OrderReWeight"]) is not list:
            print "The Reweight should be a list!"
            return False
        if para["Order"]+1 is not len(para["OrderReWeight"]):
            print "The Reweight numbers should be equal to Order!"
            return False

    def to_string(self, pid=0):
        #set Seed here so that each job has it own rng seed
        self.para["Seed"] = int(random.random()*2**30)
        return Job.to_string(self, pid)

class JobConsistLoop(Job):
    '''job subclass for self consistent loop jobs'''
    def __init__(self, para):
        Job.__init__(self, para)
        self.keep_cpu_busy = False
        self.para["Type"] = "DYSON"
        self.name = "DYSON"

    def to_string(self, pid=0):
        return Job.to_string(self, pid)

if __name__ == "__main__":
    A = JobMonteCarlo({
        "__Execute": ["python", "./monte_carlo.exe"],
        "__IsCluster":True,
        "__Duplicate":3,
        "IsForever" : True,
        "Sample" : 1000000,
        "Sweep" : 10,
        "Toss" : 1000,
        "DoesLoad" : True,
        "StartFromBare" :False,
        "Type" : 2,
        "Lx" :  4,
        "Ly" :  4,
        "Jcp" :  1.0,
        "InitialBeta" :  0.5,
        "DeltaBeta" :  0.05,
        "FinalBeta" :  0.9,
        "Order" :  2,
        "OrderReWeight" : [1,5],
        "WormSpaceReweight" : 0.5
    })
    print A.to_string(1)
    print A.execute
