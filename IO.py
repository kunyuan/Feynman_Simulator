#!/usr/bin/python
import pprint

def SaveDict(filename, mode, root):
    with open(filename, mode) as f:
        pprint.pprint(root, f)

def LoadDict(filename):
    with open(filename, "r") as f:
        content=f.read()
        return eval(content)
def Simple():
    print "Hello World!"
