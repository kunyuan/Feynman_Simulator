#!/usr/bin/python
import pprint
import cPickle as pickle
from numpy import *
#all numpy symbols have to be imported as * in order to read "array([...])" in .txt file with LoadDict function

set_printoptions(threshold=nan) #make sure numpy will print all elements, so that SaveDict and LoadDict will work even for very large array

def SaveDict(filename, mode, root):
    if filename[-4:]!=".txt":
        filename+=".txt"
    with open(filename, mode) as f:
        f.write(pprint.pformat(root))

def LoadDict(filename):
    if filename[-4:]!=".txt":
        filename+=".txt"
    with open(filename, "r") as f:
        return eval(f.read())

def SaveBigDict(filename, root):
    #TODO: better to add an file read-lock here
    if filename[-4:]!=".pkl":
        filename+=".pkl"
    with open(filename, "w") as f:
        pickle.dump(root, f, pickle.HIGHEST_PROTOCOL)

def LoadBigDict(filename):
    if filename[-4:]!=".pkl":
        filename+=".pkl"
    with open(filename, "r") as f:
        root=pickle.load(f)
    return root

def InspectDict(Dict):
    for e in Dict.keys():
        if type(Dict[e]) is ndarray:
            Dict[e]="{0} with {1}".format(type(Dict[e]), Dict[e].shape)
        elif type(Dict[e]) is dict:
            Dict[e]=InspectDict(Dict[e])
    return Dict

def InspectBigFile(filename):
    pprint.pprint(InspectDict(LoadBigDict(filename)))


