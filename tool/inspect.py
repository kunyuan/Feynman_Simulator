#!/usr/bin/env python
import sys
import numpy as np
import pprint
sys.path.append("../") #add the root dir into PYTHONPATH
import IO

def InspectDict(Dict):
    for e in Dict.keys():
        if type(Dict[e]) is dict:
            Dict[e]=InspectDict(Dict[e])
        elif type(Dict[e]) is np.ndarray:
            Dict[e]="{0} with shape {1}".format(type(Dict[e]), Dict[e].shape)
        else:
            try:
                if len(Dict[e])>100:
                   Dict[e]="{0} with shape {1}".format(type(Dict[e]), len(Dict[e]))
            except:
                pass
    return Dict

def InspectBigFile(filename):
    pprint.pprint(InspectDict(IO.LoadBigDict(filename)))

if __name__=="__main__":
    if len(sys.argv) < 2:
        sys.exit('Usage: %s file name to inspect' % sys.argv[0])
    InspectBigFile(sys.argv[1])

