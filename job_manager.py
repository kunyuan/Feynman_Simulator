#!/usr/bin/python
'''This is the job manage script.
   This is a quite universal code for all different type of simulations'''
import os, sys
import time
import subprocess
import logging
import IO
import inlist

PROCLIST = []
PROCLIST_BACK = []
workdir="."
os.system("cp IO.py "+workdir)
logging.basicConfig(filename=workdir+"/project.log",
        level=logging.INFO,
        format="\n[job.daemon][%(asctime)s][%(levelname)s]:\n%(message)s",
        datefmt='%y/%m/%d %H:%M:%S')

INFILEPATH = os.path.abspath(workdir+"/infile")
OUTFILEPATH = os.path.abspath(workdir+"/outfile")
PURE_BACK = False

class JobAtom():
    '''atom class of all jobs'''
    def __init__(self, pid, bundle):
        self.pid = pid
        execu=bundle.control["__Execute"]
        if type(execu) is str:
            self.execute = os.path.abspath(execu)
        elif type(bundle.control["__Execute"]) is list:
            self.execute = " ".join([os.path.abspath(e) if os.path.isfile(e) else e for e in execu])
        else:
            print "Jobs.execute should be a list or str!"
        self.is_cluster=bundle.control["__IsCluster"]
        self.auto_run=bundle.control["__AutoRun"]
        self.keep_cpu_busy=bundle.control["__KeepCPUBusy"]
        self.name = bundle.name
        self.para = bundle.to_dict(pid)
        return

    def get_job_name(self):
        '''get the name of JobAtom object'''
        return "Job({0}).{1}".format(self.name, self.pid)

def construct_job_queue(to_do):
    '''construct JobAtom queue from Job class '''
    logging.info("Constructing the job queue...")
    job_queue = []
    global PURE_BACK
    pid = 0
    #search folder for old jobs, the new pid=largest old pid+1
    if os.path.exists(INFILEPATH):
        filelist = [int(elem.split('.')[0].split('_')[-1]) for elem in os.listdir(INFILEPATH)]
        filelist.sort()
        if len(filelist) != 0:
            pid = filelist[-1]
    #bundle is class job
    for bundle in [e for e in to_do if e.control["__KeepCPUBusy"]== False]:
    #running the jobs doesn't use much cpu first
        for _ in range(0, bundle.control["__Duplicate"]):
            pid += 1
            if not bundle.control["__IsCluster"]:
                PURE_BACK = True
            job_queue.append(JobAtom(pid, bundle))

    for bundle in [e for e in to_do if e.control["__KeepCPUBusy"] == True]:
    #running the jobs use much cpu next
        for _ in range(0, bundle.control["__Duplicate"]):
            pid += 1
            if not bundle.control["__IsCluster"]:
                PURE_BACK = False
            job_queue.append(JobAtom(pid, bundle))

    logging.info("Constructed the job queue!")
    return job_queue

def check_status():
    ''' check the status of submitted jobs,
    if the job is done, remove it from PROCLIST so new job could be submitted'''
    for elemp in PROCLIST:
        if elemp[0].poll() is not None:
            PROCLIST.remove(elemp)
            logging.info(elemp[1].get_job_name()+" is ended!")
            print elemp[1].get_job_name()+" is ended..."
    for elemp in PROCLIST_BACK:
        if elemp[0].poll() is not None:
            PROCLIST_BACK.remove(elemp)
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
    jobname = homedir.split("/")[-1]+"."+job_atom.name

    infile = INFILEPATH+"/_in_{0}_{1}".format(job_atom.name, job_atom.pid)
    outfile = OUTFILEPATH+"/out_{0}_{1}.txt".format(
        job_atom.name, job_atom.pid)
    jobfile = os.path.abspath(workdir+"/_job_{0}_{1}.sh".format(
        job_atom.name, job_atom.pid))
    IO.SaveDict(infile, "w", job_atom.para)
    f_allinput = open(os.path.abspath(workdir+"/all_input.log"), "a")
    f_allinput.write("Job ID: {0}, Job name: {1}\n".format(
            job_atom.pid, job_atom.name))
    f_allinput.write(str(job_atom.para))
    f_allinput.close()
    if job_atom.is_cluster:
        fjob = open(jobfile, "w")
        fjob.write("#!/bin/sh\n"+"#PBS -N "+jobname+"\n")
        fjob.write("#PBS -o "+homedir+"/Output\n")
        fjob.write("#PBS -e "+homedir+"/Error\n")
        fjob.write("cd "+homedir+"\n")
        fjob.write(job_atom.execute+" "+infile)
        fjob.close()
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

def StopTheWorld():
    check_status()
    ##terminate all background process here
    for elem in PROCLIST_BACK+PROCLIST:
        elem[0].kill()
        logging.info(elem[1].get_job_name()+" is ended!")
        print elem[1].get_job_name()+" is ended..."
        sys.exit(0)

if __name__ == "__main__":
    logging.info("Jobs manage daemon is started...")
    JOBQUEUE = construct_job_queue(inlist.TO_DO)
    #print [e.keep_cpu_busy for e in JOBQUEUE]
    i = 0
    for ATOM in JOBQUEUE:
        while ATOM.is_cluster is False and len(PROCLIST)>=inlist.CPU:
            check_status()
            time.sleep(inlist.SLEEP)
        submit_job(ATOM)

    check_status()
    while len(PROCLIST+PROCLIST_BACK) != 0:
        try:
            time.sleep(inlist.SLEEP)
            check_status()
        except KeyboardInterrupt:
            StopTheWorld()

    logging.info("Jobs manage daemon is ended...")
