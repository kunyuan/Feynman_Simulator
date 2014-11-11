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
        self.para = para
        self.name = ""

    def __check_parameters__(self, para):
        if para["__Execute"] is "":
            print "Please specify the executive file name!"
            return False
        if type(para["DoesLoad"]) is not bool:
            print "DoesLoad should be a bool!"
            return False
        return True

    def key_to_string(self, key):
        '''change a key in the parameter dictionary into a string'''
        if type(self.para[key])==bool:
            if self.para[key]:
                return "1    #{0}\n".format(key)
            else:
                return "0    #{0}\n".format(key)
        elif type(self.para[key])==str:
            return self.para[key]+"    #{0}\n".format(key)
        elif type(self.para[key])==list:
            return "{0}    #{1}\n".format(",".join([str(elem)
                                 for elem in self.para[key]]), key)
        else:
            return "{0}    #{1}\n".format(self.para[key], key)

    def to_string(self, pid=0):
        '''output the corresponding string of the job class'''
        self.para["pid"] = pid
        input_str = self.key_to_string("DoesLoad")
        input_str += self.key_to_string("StartFromBare")
        input_str += self.key_to_string("pid")
        input_str += self.key_to_string("L")
        input_str += self.key_to_string("Jcp")
        input_str += self.key_to_string("initialBeta")
        input_str += self.key_to_string("deltaBeta")
        input_str += self.key_to_string("finalBeta")
        input_str += self.key_to_string("Order")
        return input_str

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
        if type(para["OrderReweight"]) is not list:
            print "The Reweight should be a list!"
            return False
        if para["Order"] is not len(para["OrderReweight"]):
            print "The Reweight numbers should be equal to Order!"
            return False

    def to_string(self, pid=0):
        input_str = self.key_to_string("Type")
        input_str = input_str+Job.to_string(self, pid)
        input_str += self.key_to_string("Toss")
        input_str += self.key_to_string("Sample")
        input_str += self.key_to_string("Sweep")
        self.para["Seed"] = -int(random.random()*2**30)
        input_str += self.key_to_string("Seed")
        input_str += self.key_to_string("WormSpaceReweight")
        input_str += self.key_to_string("OrderReweight")
        return input_str

class JobConsistLoop(Job):
    '''job subclass for self consistent loop jobs'''
    def __init__(self, para):
        Job.__init__(self, para)
        self.keep_cpu_busy = False
        self.para["Type"] = "DYSON"
        self.name = "DYSON"

    def to_string(self, pid=0):
        input_str = self.key_to_string("Type")
        input_str = input_str+Job.to_string(self, pid)
        input_str += self.key_to_string("OrderAccepted")
        input_str += self.key_to_string("ErrorThreshold")
        input_str += self.key_to_string("SleepTime")
        return input_str

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
        "iniBeta" :  0.5,
        "dBeta" :  0.05,
        "finalBeta" :  0.9,
        "Order" :  2,
        "Reweight" : [1,5],
        #"ReadFile" : "0.90_1_coll",
        "WormSpaceReweight" : 0.5
    })
    print A.to_string(1)
    print A.execute
