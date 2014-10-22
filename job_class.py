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
        if type(para["IsLoad"]) is not bool:
            print "IsLoad should be a bool!"
            return False
        return True

    def key_to_string(self, key):
        '''change a key in the parameter dictionary into a string'''
        if type(self.para[key])==bool:
            if self.para[key]:
                return ".true.    #{0}\n".format(key)
            else:
                return ".false.    #{0}\n".format(key)
        elif type(self.para[key])==str:
            return self.para[key]+"\n"
        elif type(self.para[key])==list:
            return "{0}    #{1}\n".format(",".join([str(elem)
                                 for elem in self.para[key]]), key)
        else:
            return "{0}    #{1}\n".format(self.para[key], key)

    def to_string(self, pid=0):
        '''output the corresponding string of the job class'''
        self.para["pid"] = pid
        input_str = self.key_to_string("pid")
        input_str += self.key_to_string("L")
        input_str += self.key_to_string("Jcp")
        input_str += self.key_to_string("iniBeta")
        input_str += self.key_to_string("dBeta")
        input_str += self.key_to_string("finalBeta")
        input_str += self.key_to_string("Order")
        input_str += self.key_to_string("IsLoad")
        return input_str

class JobMonteCarlo(Job):
    '''job subclass for monte carlo jobs'''
    def __init__(self, para):
        Job.__init__(self, para)
        self.keep_cpu_busy = True
        self.para["Type: MC"] = 2
        self.name = "MC"

    def __check_parameters__(self, para):
        if Job.__check_parameters__(self, para) is False:
            return False
        if type(para["Reweight"]) is not list:
            print "The Reweight should be a list!"
            return False
        if para["Order"] is not len(para["Reweight"]):
            print "The Reweight numbers should be equal to Order!"
            return False

    def to_string(self, pid=0):
        input_str = Job.to_string(self, pid)
        input_str += self.key_to_string("Type: MC")
        input_str += self.key_to_string("Toss")
        input_str += self.key_to_string("Sample")
        input_str += self.key_to_string("Sweep")
        self.para["Seed"] = -int(random.random()*2**30)
        input_str += self.key_to_string("Seed")
        input_str += self.key_to_string("ReadFile")
        input_str += self.key_to_string("Worm/Norm")
        input_str += self.key_to_string("Reweight")
        return input_str

class JobConsistLoop(Job):
    '''job subclass for self consistent loop jobs'''
    def __init__(self, para):
        Job.__init__(self, para)
        self.keep_cpu_busy = False
        self.para["Type: SCL"] = 1
        self.name = "SCL"

    def to_string(self, pid=0):
        input_str = Job.to_string(self, pid)
        input_str += self.key_to_string("Type: SCL")
        input_str += self.key_to_string("ReadFile")
        return input_str

class JobIntegration(Job):
    '''job subclass for numerical integration jobs'''
    def __init__(self, para):
        Job.__init__(self, para)
        self.keep_cpu_busy = True
        self.para["Type: NI"] = 3
        self.name = "NI"

    def to_string(self, pid=0):
        input_str = Job.to_string(self, pid)
        input_str += self.key_to_string("Type: NI")
        return input_str

class JobOutputOrder(Job):
    '''job subclass for output different order Sigma and Chi'''
    def __init__(self, para):
        Job.__init__(self, para)
        self.keep_cpu_busy = True
        self.para["Type: OO"] = 5
        self.name = "OO"

    def to_string(self, pid=0):
        input_str = Job.to_string(self, pid)
        input_str += self.key_to_string("Type: OO")
        input_str += self.key_to_string("ReadFile")
        return input_str

class JobOutputLoop(Job):
    '''job subclass for output loop jobs'''
    def __init__(self, para):
        Job.__init__(self, para)
        self.keep_cpu_busy = False
        self.para["Type: OL"] = 4
        self.name = "OL"

    def to_string(self, pid=0):
        input_str = Job.to_string(self, pid)
        input_str += self.key_to_string("Type: OL")
        input_str += self.key_to_string("ReadFile")
        return input_str

class JobDebug(Job):
    '''job subclass for debug jobs'''
    def __init__(self, para):
        Job.__init__(self, para)
        self.keep_cpu_busy = False
        self.para["Type: BG"] = 6
        self.name = "BG"

    def to_string(self, pid=0):
        input_str = Job.to_string(self, pid)
        input_str += self.key_to_string("Type: BG")
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
        "IsLoad" : True,
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
        "Worm/Norm" : 0.5
    })
    print A.to_string(1)
    print A.execute
