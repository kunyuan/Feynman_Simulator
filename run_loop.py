#!/usr/bin/python
'''This code could be used to run self consistent loop or output loop'''
import os
import sys
import subprocess
import time
import logging

INTERVAL =300 
EXEC = "./src/FeynmanSimulator.exe"

def run_loop(infile):
    '''the loop to do self consisent calculation or output variables'''
    homedir = os.getcwd()
    logging.basicConfig(filename=homedir+"/data/project.log", level=logging.INFO,
         format="\n[loop.daemon][%(asctime)s][%(levelname)s]:\n%(message)s",
         datefmt='%y/%m/%d %H:%M:%S')
    execf = os.path.abspath(EXEC)
    logging.info("Loop daemon started!")
    #logging.info(title+" is the target!")
    print "loop daemon started..."
    #time.sleep(INTERVAL)
    i = 0
    while True:
        i += 1
        logging.info("Loop "+str(i)+" running...")
        os.system("rm read_list.dat")
        os.system("cat readfile/*.dat >>read_list.dat")
        try:
            proc = subprocess.Popen(homedir+"/collapse_data.exe")
            exitcode = proc.wait()
            if exitcode is not 0:
                raise RuntimeError('collapse_data return a value '+exitcode)
        except subprocess.CalledProcessError as err:
            ret = err.returncode
            logging.error('collapse_data.exe return a non-zero value '
                +str(ret)+', something happened!')
        except RuntimeError as err:
            logging.error(err.message)
        except:
            logging.error('Collapse_data failed!')
        else:
            try:
                #the execf process share the same stdout as the run_loop.py has
                os.system('exec '+execf+' '+infile)
            except:
                logging.error('loop error!')

        logging.info("Loop "+str(i)+" done!")
        logging.info("Sleeping...")
        time.sleep(INTERVAL)
        logging.info("Loop daemon ended!")

INFILE = sys.argv[1]
INFILE = os.path.abspath(INFILE)

print "Start loop..."
run_loop(INFILE)
