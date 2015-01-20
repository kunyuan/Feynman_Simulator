#!/usr/bin/env python
from logger import *
import sys
import pprint

def Load(FileName):
    log.info("Loading Parameters...")
    d=IO.LoadDict(FileName)
    log.info("Loaded parameters:\n"+pprint.pformat(d))
    Assert(d["Job"]["Type"]=="DYSON", "The job type should be DYSON, not {0}".format(d["Job"]["Type"]))
    return d["Job"], d["Para"]

def Save(para, FileName, Mode="a"):
    log.info("Saving Parameters...")
    root={"Para":para}
    IO.SaveDict(FileName, Mode, root)

def GetVersion(MessageFile):
    try:
        Version=IO.LoadDict(MessageFile)["Version"]
    except:
        Version=0
    finally:
        return Version

def BroadcastMessage(MessageFile, Dict):
    log.info("Broadcast Message")
    IO.SaveDict(MessageFile, "w", Dict)

if __name__=="__main__":
    p=Load("../infile/_in_DYSON_1")
    Save(p, "test","w")
    print p["Dyson"]["Order"]



