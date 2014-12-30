#!/usr/bin/python
import pprint
import hickle
from numpy import *
set_printoptions(threshold=nan) #make sure numpy will print all elements, so that SaveDict and LoadDict will work even for very large array

def SaveDict(filename, mode, keystr, root):
    with open(filename, mode) as f:
        f.write(keystr+"="+pprint.pformat(root))

def LoadDict(filename, keystr):
    with open(filename, "r") as f:
        content=f.read()
        exec(content)
        return locals()[keystr]

def SaveBigDict(filename, root):
    with open(filename+".hkl", "w") as f:
        hickle.dump(root, f)

def LoadBigDict(filename):
    with open(filename+".hkl", "w") as f:
        root=hickle.load(f)
    return root




