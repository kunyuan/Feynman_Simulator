#!/usr/bin/env python
from logger import *
import sys
import pprint

def Load(FileName):
    log.info("Loading Parameters...")
    d=IO.LoadDict(FileName)
    #log.info("Loaded parameters:\n"+pprint.pformat(d))
    Assert(d["Job"]["Type"]=="DYSON", "The job type should be DYSON, not {0}".format(d["Job"]["Type"]))
    return d["Job"], d["Para"]

def LoadPara(FileName):
    log.info("Loading Parameters...")
    d=IO.LoadDict(FileName)
    #log.info("Loaded parameters:\n"+pprint.pformat(d))
    return d["Para"]

def Save(FileName, para):
    log.info("Saving Parameters...")
    root={"Para":para}
    IO.SaveDict(FileName, "w", root)

def BroadcastMessage(MessageFile, Dict):
    log.info("Broadcast Message")
    IO.SaveDict(MessageFile, "w", Dict)

if __name__=="__main__":
    p=Load("../infile/_in_DYSON_1")
    Save(p, "test","w")
    print p["Dyson"]["Order"]



