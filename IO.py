#!/usr/bin/python
import pprint
import cPickle as pickle
from numpy import *

set_printoptions(threshold=nan) #make sure numpy will print all elements, so that SaveDict and LoadDict will work even for very large array

def SaveDict(filename, mode, root):
    with open(filename+".txt", mode) as f:
        f.write(pprint.pformat(root))

def LoadDict(filename):
    with open(filename+".txt", "r") as f:
        return eval(f.read())

def SaveBigDict(filename, root):
    with open(filename+".pkl", "w") as f:
        pickle.dump(root, f, pickle.HIGHEST_PROTOCOL)

def LoadBigDict(filename):
    with open(filename+".pkl", "r") as f:
        root=pickle.load(f)
    return root

