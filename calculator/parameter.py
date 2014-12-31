#!/usr/bin/env python
from logger import *
import sys
import pprint
sys.path.append("../") #add the root dir into PYTHONPATH
import IO

class Parameter:
    def Load(self, FileName):
        log.info("Loading Parameters...")
        d=IO.LoadDict(FileName)
        self.__dict__=d["Para"]
        print self.__dict__
        log.info("Loaded parameters:\n"+pprint.pformat(self.__dict__))

    def Save(self, FileName, Mode="a"):
        log.info("Saving Parameters...")
        with open(FileName, Mode) as f:
            f.write("Para="+pprint.pformat(self.__dict__))

if __name__=="__main__":
    p=Parameter()
    p.Load("../data/infile/_in_DYSON_1")
    p.Save("test.txt","w")
    print p.Order



