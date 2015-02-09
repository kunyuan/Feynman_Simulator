#!/usr/bin/env python
import os
filename="id_job.log"
lines = [line.rstrip('\n') for line in open(filename)]
for jobid in lines:
    print "Kill {0}".format(jobid)
    os.system("qdel {0}".format(jobid))
    #os.system("qsig -s SIGINT {0}".format(jobid))
