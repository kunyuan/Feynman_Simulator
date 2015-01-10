#!/usr/bin/env python
import numpy as np
from logger import *
import os, sys, model, weight

StatisFilePattern="_statistics"
FileList = [f for f in os.listdir(workspace) if os.path.isfile(os.path.join(workspace,f))]
StatisFileList=[f for f in FileList if f.find(StatisFilePattern) is not -1]

if __name__=="__main__":
    print StatisFileList

