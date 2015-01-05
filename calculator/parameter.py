#!/usr/bin/env python
from logger import *
import sys
import pprint
sys.path.append("../") #add the root dir into PYTHONPATH
import IO

def Load(FileName):
    log.info("Loading Parameters...")
    d=IO.LoadDict(FileName)["Para"]
    log.info("Loaded parameters:\n"+pprint.pformat(d))
    Assert(d["Job"]["Type"]=="DYSON", "The job type should be DYSON, not {0}".format(d["Job"]["Type"]))
    return d

def Save(para, FileName, Mode="a"):
    log.info("Saving Parameters...")
    root={"Para":para}
    IO.SaveDict(FileName, Mode, root)

if __name__=="__main__":
    p=Load("../data/infile/_in_DYSON_1")
    Save(p, "test","w")
    print p["Dyson"]["Order"]



