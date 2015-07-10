#!/usr/bin/env python
import os
import argparse
from multiprocessing import Pool

def delete(jobid):
    print "Deleting job {0}".format(jobid)
    os.system("qdel {0}".format(jobid))
def terminate(jobid):
    print "Terminating job {0}".format(jobid)
    os.system("qsig -s SIGTERM {0}".format(jobid))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--force", action='store_true', help="kill jobs with qdel command")
    parser.add_argument('-r', '--range', nargs=2, required=False, type=int, help="--range start end: kill [start, end) jobs")
    args = parser.parse_args()

    if args.range is not None:
        lines=range(args.range[0],args.range[1]) 
    else:
        filename="id_job.log"
        lines = [line.rstrip('\n') for line in open(filename)]
    for i in lines:
        if args.force:
            delete(i) 
        else:
            terminate(i)
