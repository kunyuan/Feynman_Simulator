#!/usr/bin/env python
from logger import *

class Parameter:
    def __formator__(self,para):
        for (k,v) in para.items():
            if type(v) is str:
                para[k]="'"+v+"'"
        return "\n".join([k+" = "+str(v) for (k,v) in para.items()])

    def Load(self, FileName):
        log.info("Loading Parameters...")
        with open(FileName, "r") as f:
            lines=f.read().splitlines()
        for e in lines:
            exec("self."+e)
        log.info("Loaded parameters:\n"+self.__formator__(self.__dict__))

    def Save(self, FileName, Mode="a"):
        log.info("Saving Parameters...")
        with open(FileName, Mode) as f:
            f.write(self.__formator__(self.__dict__))

if __name__=="__main__":
    p=Parameter()
    p.Load("../data/infile/_in_DYSON_1")
    p.Save("test.txt","w")
    print p.Order



