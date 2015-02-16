#!/usr/bin/env python
import os
import argparse
filename="id_job.log"
lines = [line.rstrip('\n') for line in open(filename)]
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--force", help="use file path to find the input file")
args = parser.parse_args()

for jobid in lines:
    print "Killing {0}".format(jobid)
    if args.force:
        os.system("qdel {0}".format(jobid))
    else:
        os.system("qsig -s SIGTERM {0}".format(jobid))
