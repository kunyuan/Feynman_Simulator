#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Last modified: 
import sys
import numpy as np
import parameter as para
import weight, plot
import parameter
from logger import *


def Vec2Index(vec, _map): 
  Index = vec[0]
  for i in range(1,_map.Dim):
      Index = Index * _map.L[i] + vec[i]
  return Index 


def Index2Vec(index, _map):
  v = np.zeros(_map.Dim)
  for i in reversed(range(1, _map.Dim)):
      v[i] = index % _map.L[i]
      index = index / _map.L[i]
  v[0] = index
  return v

def Shift(vec, _map):
    for i in range(_map.Dim):
        if vec[i] < 0:
            vec[i] += _map.L[i]
        if vec[i] >= _map.L[i]:
            vec[i] -= _map.L[i]
    return vec
        
def CoordiIndex(index_out, index_in, _map):
    vec_out = Index2Vec(index_out, _map)
    vec_in = Index2Vec(index_in, _map)
    vec = vec_out - vec_in
    vec = Shift(vec, _map)
    return Vec2Index(vec, _map)


if __name__ == "__main__":
    global ParaFile
    ParaFile="0_DYSON_para" 
    para=parameter.LoadPara(ParaFile)

    WeightPara={"NSublat": para["Lattice"]["NSublat"], "L":para["Lattice"]["L"],
                "Beta": float(para["Tau"]["Beta"]), "MaxTauBin": para["Tau"]["MaxTauBin"]}
    Map=weight.IndexMap(**WeightPara)

    index_out = 34
    index_in = 56

    print CoordiIndex(index_out, index_in, Map)
    
