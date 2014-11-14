#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import os, sys
lib_path = os.path.abspath('../data')
sys.path.append(lib_path)
import lattice

color=('r','g','b')
points=lattice.points
for coord, label, sub in points:
    x,y=coord;
    plt.scatter(x,y,s=100,c=color[sub])
    plt.annotate(
            str(label),
            xy = (x, y), xytext = (15, 10),
            textcoords = 'offset points', ha = 'right', va = 'bottom',
            arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))

import interaction
lines=interaction.line
for start,end, sub in lines:
    x,y=zip(start,end)
    plt.plot(x,y,c=color[sub],lw=2)

plt.show()
