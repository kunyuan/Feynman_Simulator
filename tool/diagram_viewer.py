#!/usr/bin/env pythonw
#"key_press_event" in MacOSX requires pythonw!!!
import os
import re
import subprocess
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import matplotlib.image as mpimg

path="../data/diagram"
figformat="jpg"
engine=["dot","neato","sfdp","fdp"]

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

def GetImage(filename, graph_engine=engine[0]):
    imagefile=filename[:-2]+figformat
    shellstr=graph_engine+" -Tjpg -Gsize=10,15\! -Gdpi=80 "+filename+" -o "+imagefile
    os.system(shellstr)
    return imagefile
    
walk=Enumerator()
fname=walk(0)
f,ax = plt.subplots(1,2,figsize=(5,5))
ax[0].axis('off')
ax[1].axis('off')
img=mpimg.imread(GetImage(fname,engine[0]))
imgplot = ax[0].imshow(img)
text=ax[1].text(0.5, 0.5,GetComment(fname),
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax[1].transAxes)
title=f.text(0.5, 0.98,fname+"\nengine:"+engine[0],
     horizontalalignment='center',
     verticalalignment='center',
     fontsize=14, fontweight='bold', alpha=0.5)

class Index:
    def __init__(self):
        self._engine=0

    def draw(self, fname):
        img=mpimg.imread(GetImage(fname,engine[self._engine]))
        imgplot.set_data(img)
        imgplot.set_extent((0,img.shape[1],0,img.shape[0]))
        text.set_text(GetComment(fname))
        title.set_text(fname+"\nengine:"+engine[self._engine])
        plt.draw()

    def next(self, event):
        self.draw(walk(1))

    def prev(self, event):
        self.draw(walk(-1))
    
    def retry(self, event):
        self._engine += 1
        if self._engine>=len(engine):
                self._engine=0
        self.draw(walk(0))
    
    def key(self, event):
        if event.key=='j' or event.key=='l' or event.key=='n' or event.key=='right' or event.key=='down':
            self.next(event)
        if event.key=='k' or event.key=='h' or event.key=='p' or event.key=='left' or event.key=='up':
            self.prev(event)
        if event.key=='g':
            self.retry(event)

callback = Index()
axretry = plt.axes([0.59, 0.01, 0.1, 0.05])
axprev = plt.axes([0.7, 0.01, 0.1, 0.05])
axnext = plt.axes([0.81, 0.01, 0.1, 0.05])
bretry=Button(axretry,'Feel lucky')
bretry.on_clicked(callback.retry)
bnext = Button(axnext, 'Next')
bnext.on_clicked(callback.next)
bprev = Button(axprev, 'Previous')
bprev.on_clicked(callback.prev)
f.canvas.mpl_connect("key_press_event",callback.key)

f.canvas.set_window_title("'j,l,n,right,down' to next, 'k,h,p,left,up' to prev, 'g' to change engine") 
plt.tight_layout(pad=0.0, w_pad=0.0, h_pad=0.0)
plt.show()
