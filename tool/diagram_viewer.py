#!/usr/bin/env python
import os
import re
import subprocess
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import matplotlib.image as mpimg

path="../data/diagram"
f,ax = plt.subplots(1,2,figsize=(5,5))
ax[0].axis('off')
ax[1].axis('off')

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

class Enumerator:
    def list_files(self):
        self.files = []
        for name in [e for e in os.listdir(path) if e[-3:]==".gv"]:
            self.files.append(os.path.join(path, name))
        self.files=natural_sort(self.files)

    def __init__(self):
        self.list_files()
        self.counter = 0

    def __call__(self,num):
        self.counter = self.counter + num
        if self.counter>=len(self.files):
            self.counter=0
        elif self.counter<0:
            self.counter=len(self.files)-1
        return self.files[self.counter]

def GetComment(filename):
    with open(filename,"r") as f:
        flist=[]
        for line in f:
            if line[0:2]=="//":
                flist.append(line[2:])
    return "\n".join(flist)

def GetImage(filename):
    imagefile=filename[:-2]+"png"
    shellstr="dot -Tpng "+filename+" -o "+imagefile
    os.system(shellstr)
    return imagefile
    
walk=Enumerator()
fname=walk(0)
img=mpimg.imread(GetImage(fname))
imgplot = ax[0].imshow(img)
text=ax[1].text(0.5, 0.5,GetComment(fname),
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax[1].transAxes)
title=plt.gcf().text(0.5, 0.98,fname,
     horizontalalignment='center',
     verticalalignment='center',
     fontsize=14, fontweight='bold', alpha=0.5)

class Index:
    ind = 0
    def next(self, event):
        fname=walk(1)
        img=mpimg.imread(GetImage(fname))
        imgplot.set_data(img)
        text.set_text(GetComment(fname))
        title.set_text(fname)
        plt.draw()

    def prev(self, event):
        fname=walk(-1)
        img=mpimg.imread(GetImage(fname))
        imgplot.set_data(img)
        text.set_text(GetComment(fname))
        title.set_text(fname)
        plt.draw()

callback = Index()
axprev = plt.axes([0.7, 0.05, 0.1, 0.075])
axnext = plt.axes([0.81, 0.05, 0.1, 0.075])
bnext = Button(axnext, 'Next')
bnext.on_clicked(callback.next)
bprev = Button(axprev, 'Previous')
bprev.on_clicked(callback.prev)

plt.tight_layout(pad=0.0, w_pad=0.0, h_pad=0.0)
plt.show()
