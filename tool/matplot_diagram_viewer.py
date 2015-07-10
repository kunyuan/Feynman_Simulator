#!/usr/bin/env python
import os
import re
import subprocess
import matplotlib
matplotlib.use('Qt4Agg')
#"key_press_event" does not work for MacOSX backend!!!
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import matplotlib.image as mpimg

path="../data/diagram"
#"jpg" or "png"
figformat="png"
engine=["dot","neato","sfdp","fdp"]
next_keys=('j','l','n','right','down')
prev_keys=('k','h','p','left','up')
retry_keys=('c','g')

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

class Enumerator:
    def __init__(self):
        self.counter = 0
        self.files=natural_sort([e[:-3] for e in os.listdir(path) if e[-3:]==".gv"])

    def __call__(self,num):
        self.counter = self.counter + num
        if self.counter>=len(self.files):
            self.counter=0
        elif self.counter<0:
            self.counter=len(self.files)-1
        return os.path.join(path, self.files[self.counter])
            
    def CurrentFile(self):
        return self.files[self.counter]

def GetComment(filename):
    '''filename has no .gv'''
    with open(filename+".gv","r") as f:
        flist=[]
        for line in f:
            if line[0:2]=="//":
                flist.append(line[2:])
    return "\n".join(flist)

def GetImage(filename, graph_engine=engine[0]):
    '''filename has no .gv'''
    imagefile=filename+"."+figformat
    print imagefile
    shellstr=graph_engine+" -T"+figformat+" -Gsize=10,15\! -Gdpi=80 "+filename+".gv -o "+imagefile
    os.system(shellstr)
    return imagefile
    
walk=Enumerator()
fname=walk(0)
f,ax = plt.subplots(1,2,figsize=(5,5))
ax[0].axis('off')
ax[1].axis('off')
img=mpimg.imread(GetImage(fname,engine[0]))
f.canvas.set_window_title(fname+".gv (engine: "+engine[0]+")")
imgplot = ax[0].imshow(img)
text=ax[1].text(0.5, 0.99,GetComment(fname),
     horizontalalignment='center',
     verticalalignment='top',
     transform = ax[1].transAxes)
status=f.text(0.4, 0.98, "", fontweight='bold', fontsize='14', alpha=0.5)
help="next:"+str(next_keys)+"; prev:"+str(prev_keys)+"; feel lucky:"+str(retry_keys)
f.text(0.05, 0.01, help, horizontalalignment='left', fontsize=11, fontweight='bold',alpha=0.5)

class Index:
    def __init__(self):
        self._engine=0

    def draw(self, fname):
        img=mpimg.imread(GetImage(fname,engine[self._engine]))
        imgplot.set_data(img)
        imgplot.set_extent((0,img.shape[1],0,img.shape[0]))
        text.set_text(GetComment(fname))
        f.canvas.set_window_title(fname+".gv (engine: "+engine[self._engine]+")")
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
    
    def save(self, event):
        fname=walk.CurrentFile()+"_"+engine[self._engine]+".pdf"
        status.set_text("Save "+fname)
        plt.savefig(fname)
        status.set_text("")
    
    def key(self, event):
        if event.key in next_keys:
            self.next(event)
        if event.key in prev_keys:
            self.prev(event)
        if event.key in retry_keys:
            self.retry(event)

callback = Index()
axsave = plt.axes([0.53, 0.01, 0.1, 0.05])
axretry = plt.axes([0.64, 0.01, 0.1, 0.05])
axprev = plt.axes([0.75, 0.01, 0.1, 0.05])
axnext = plt.axes([0.86, 0.01, 0.1, 0.05])
bsave=Button(axsave,'Save')
bsave.on_clicked(callback.save)
bretry=Button(axretry,'Feel lucky')
bretry.on_clicked(callback.retry)
bnext = Button(axnext, 'Next')
bnext.on_clicked(callback.next)
bprev = Button(axprev, 'Previous')
bprev.on_clicked(callback.prev)
f.canvas.mpl_connect("key_press_event",callback.key)

plt.tight_layout(pad=0.0, w_pad=0.0, h_pad=0.0)
plt.show()
