#!/usr/bin/env python
'''This is the job manage script.
   This is a quite universal code for all different type of simulations'''
import os, sys, copy, signal
import time
import subprocess
import logging
import IO
import inlist

PROCLIST = []
PROCLIST_BACK = []
workdir=os.path.abspath(".")
logging.basicConfig(filename=workdir+"/project.log",
        level=logging.INFO,
        format="\n[job.daemon][%(asctime)s][%(levelname)s]:\n%(message)s",
        datefmt='%y/%m/%d %H:%M:%S')

INFILEPATH = os.path.join(workdir,"infile")
OUTFILEPATH = os.path.join(workdir,"outfile")

class JobAtom():
    '''atom class of all jobs'''
    def __init__(self, control, pid, para):
        self.pid = pid
        self.name = control["__Type"]
        execu=control["__Execute"]
        if type(execu) is str:
            self.execute = os.path.abspath(execu)
        elif type(control["__Execute"]) is list:
            self.execute = " ".join([os.path.abspath(e) if os.path.isfile(e) else e for e in execu])
        else:
            print "Jobs.execute should be a list or str!"
        self.is_cluster=control["__IsCluster"]
        self.auto_run=control["__AutoRun"]
        self.keep_cpu_busy=control["__KeepCPUBusy"]
        if control.has_key("__PBSCommand"):
            self.pbs_command=control["__PBSCommand"]
        self.para = para
        return

    def get_job_name(self):
        '''get the name of JobAtom object'''
        return "Job({0}).{1}".format(self.name, self.pid)

def construct_job_queue(to_do):
    '''construct JobAtom queue from Job class '''
    logging.info("Constructing the job queue...")
    job_queue = []
    job_queue_back = []
    for JobClass in to_do:
        try:
            while True:
                pid, para=JobClass.to_dict()
                atom=JobAtom(JobClass.control, pid, para)
                if atom.keep_cpu_busy:
                    job_queue.append(atom)
                else:
                    job_queue_back.append(atom)
        except IndexError:
            pass 
    logging.info("Constructed the job queue!")
    return job_queue, job_queue_back

def check_status():
    ''' check the status of submitted jobs,
    if the job is done, remove it from PROCLIST so new job could be submitted'''
    PROC=(PROCLIST, PROCLIST_BACK)
    #print PROC
    for plist in PROC:
        for elemp in plist:
            if elemp[0].poll() is not None:
                plist.remove(elemp)
                logging.info(elemp[1].get_job_name()+" is ended!")
                print elemp[1].get_job_name()+" is ended..."
    return

def submit_job(job_atom):
    '''submit a job to cluster or your own computer'''

    if os.path.exists(INFILEPATH) is not True:
        os.system("mkdir "+INFILEPATH)
    if os.path.exists(OUTFILEPATH) is not True:
        os.system("mkdir "+OUTFILEPATH)

    homedir = os.getcwd()
    _, tail = os.path.split(homedir)
    jobname = tail+"."+job_atom.name

    infile = os.path.join(INFILEPATH,"_in_{0}_{1}".format(job_atom.name, job_atom.pid))
    outfile = os.path.join(OUTFILEPATH,"out_{0}_{1}.txt".format(job_atom.name, job_atom.pid))
    jobfile = os.path.join(workdir,"_job_{0}_{1}.sh".format(job_atom.name, job_atom.pid))
    IO.SaveDict(infile, "w", job_atom.para)
    if job_atom.is_cluster:
        with open(jobfile, "w") as fjob:
            fjob.write("#!/bin/sh\n"+"#PBS -N "+jobname+"\n")
            if hasattr(job_atom, "pbs_command"):
                fjob.write(job_atom.pbs_command+"\n")
            fjob.write("#PBS -o "+homedir+"/Output\n")
            fjob.write("#PBS -e "+homedir+"/Error\n")
            fjob.write("echo $PBS_JOBID >>"+homedir+"/id_job.log\n")
            fjob.write("cd "+homedir+"\n")
            fjob.write(job_atom.execute+" -f "+infile)
        if job_atom.auto_run:
            os.system("qsub "+jobfile)
            os.system("rm "+jobfile)
            logging.info(job_atom.get_job_name()+" submitted!")
        else:
            print "You have to run "+job_atom.get_job_name()+" by yourself!"
    else:
        if job_atom.auto_run:
            #shellstr = "exec "+job_atom.execute+" -f "+infile+" >> "+outfile
            shellstr = "exec "+job_atom.execute+" -f "+infile
            proc = subprocess.Popen(shellstr, shell=True)
            if job_atom.keep_cpu_busy:
                PROCLIST.append((proc, job_atom))
            else:
                PROCLIST_BACK.append((proc, job_atom))

            logging.info(job_atom.get_job_name()+" is started...")
            logging.info("input:\n"+str(job_atom.para))
            logging.info("PID:{0}\n".format(proc.pid))
            print job_atom.get_job_name()+" is started..."
        else:
            print "You have to run "+job_atom.get_job_name()+" by yourself!"
    return

def Exit():
    logging.info("Jobs manage daemon is ended")
    print("Jobs manage daemon is ended")
    sys.exit(0)

def StopTheWorld(signum, frame):
    """Stop the sub-process *child* if *signum* is SIGTERM. Then terminate."""
    check_status()
    ##terminate all background process here
    try:
        for elem in PROCLIST_BACK+PROCLIST:
            elem[0].terminate()
            elem[0].wait()
            logging.info(elem[1].get_job_name()+" is ended!")
            print elem[1].get_job_name()+" is ended!"
    except:
        traceback.print_exc()

if __name__ == "__main__":
    logging.info("Jobs manage daemon is started...")
    try:
        JOBQUEUE, JOBQUEUE_BACK = construct_job_queue(inlist.TO_DO)
        #print [e.keep_cpu_busy for e in JOBQUEUE]
        i = 0
        signal.signal(signal.SIGINT, StopTheWorld)
        signal.signal(signal.SIGTERM, StopTheWorld)
        for ATOM in JOBQUEUE_BACK:
            submit_job(ATOM)
        for ATOM in JOBQUEUE:
            while ATOM.is_cluster is False and len(PROCLIST)>=inlist.CPU:
                check_status()
                time.sleep(inlist.SLEEP)
            submit_job(ATOM)
        while len(PROCLIST+PROCLIST_BACK) != 0:
            check_status()
            time.sleep(inlist.SLEEP)
    except:
        traceback.print_exc()
    finally:
        Exit()
