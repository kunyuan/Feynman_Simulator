#!/usr/bin/env python
import logging
import sys
workspace="../"
sys.path.append(workspace) #add the root dir into PYTHONPATH
import IO

log = logging.getLogger()
log.setLevel(logging.INFO)
ch = logging.StreamHandler(sys.stdout)
fh = logging.FileHandler(workspace+'project.log')
ch.setLevel(logging.INFO)
fh.setLevel(logging.INFO)
formatter = logging.Formatter(fmt="[calc][%(asctime)s][%(levelname)s]:\n%(message)s",
        datefmt='%y/%m/%d %H:%M:%S')
ch.setFormatter(formatter)
fh.setFormatter(formatter)
log.addHandler(ch)
log.addHandler(fh)

def Assert(condition, info):
    if not condition:
        log.error(info)
        raise AssertionError
