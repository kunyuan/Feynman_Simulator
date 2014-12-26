#!/usr/bin/python
import pprint
import h5py
import numpy

def SaveDict(filename, mode, key, root):
    with open(filename, mode) as f:
        f.write(key+"="+pprint.pformat(root))

def LoadDict(filename, keystr):
    with open(filename, "r") as f:
        content=f.read()
        exec(content)
        return locals()[keystr]

def SaveBigDict(filename, mode, key, root):


