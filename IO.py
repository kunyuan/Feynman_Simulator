#!/usr/bin/python
import pprint

def SaveDict(filename, mode, key, root):
    with open(filename, mode) as f:
        f.write(key+"="+pprint.pformat(root))

def LoadDict(filename, keystr):
    with open(filename, "r") as f:
        content=f.read()
        exec(content)
        return locals()[keystr]
def Simple():
    print "Hello World!"
