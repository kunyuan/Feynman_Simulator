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
        self.control=para.pop("Control")
        self.pid = 0
        self.name = ""
        self.para = para

    def to_dict(self, pid=0):
        '''output the corresponding string of the job class'''
        self.para["Job"]["PID"] = pid
        self.para["Job"]["WeightFile"]="Weight"
        self.para["Job"]["MessageFile"]="Message"
        self.__set_model_specific__()
        para_={}
        para_["Para"]=self.para
        return para_

    def __check_parameters__(self, para):
        if para["Control"]["__Execute"] is "":
            print "Please specify the executive file name!"
            return False

    def __set_model_specific__(self):
        PI=3.141592653589793238
        if self.para["Model"]["Name"]=="J1J2":
            self.para["Model"]["Hopping"]=[0.0,]
            mu=1j*PI/2.0/self.para["Tau"]["Beta"]
            self.para["Model"]["ChemicalPotential"]=[mu,mu]

class JobMonteCarlo(Job):
    '''job subclass for monte carlo jobs'''
    def __init__(self, para):
        Job.__init__(self, para)
        self.para["Job"]["Type"] = "MC"
        self.name = "MC"

    def __check_parameters__(self, para):
        if Job.__check_parameters__(self, para) is False:
            return False
        if type(para["Markov"]["OrderReWeight"]) is not list:
            print "The Reweight should be a list!"
            return False
        if para["Markov"]["Order"]+1 is not len(para["Markov"]["OrderReWeight"]):
            print "The Reweight numbers should be equal to Order!"
            return False

    def to_dict(self, pid=0):
        #set Seed here so that each job has it own rng seed
        self.para["Markov"]["Seed"] = int(random.random()*2**30)
        return Job.to_dict(self, pid)

class JobDyson(Job):
    '''job subclass for self consistent loop jobs'''
    def __init__(self, para):
        Job.__init__(self, para)
        self.para["Job"]["Type"] = "DYSON"
        self.name = "DYSON"

    def to_dict(self, pid=0):
        return Job.to_dict(self, pid)

