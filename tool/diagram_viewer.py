#!/usr/bin/env python
import sys
import os
import re
from PyQt4 import QtGui,QtCore,QtSvg

diagram_input_path="../diagram/"
diagram_export_path="./diagram/"
icon="./icon/"
engine=["dot","neato","sfdp","fdp"]
qt=QtCore.Qt
prev_keys={'p': qt.Key_P,'right': qt.Key_Left, 'h':qt.Key_H}
next_keys={'n':  qt.Key_N,'left':  qt.Key_Right, 'l': qt.Key_L}
retry_keys={'space': qt.Key_Space,'g': qt.Key_G, 'f':qt.Key_F}
save_keys={'s': qt.Key_S}
zoomin_keys={'up': qt.Key_Up}
zoomout_keys={'down': qt.Key_Down}
reset_keys={'r': qt.Key_R}

class Enumerator:
    def __init__(self):
        self.counter = 0
        convert = lambda text: int(text) if text.isdigit() else text.lower() 
        alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
        self.files=sorted([e[:-3] for e in os.listdir(diagram_input_path) if e[-3:]==".gv"], \
                          key=alphanum_key)

    def __call__(self,num, IsFullPath=True):
        self.counter += num
        if self.counter>=len(self.files):
            self.counter=0
        elif self.counter<0:
            self.counter=len(self.files)-1
        if IsFullPath is True:
            return os.path.join(diagram_input_path, self.files[self.counter])
        else:
            return self.files[self.counter]

    def jumpto(self, fname, IsFullPath=True):
        self.counter=self.files.index(fname)
        return self(0)

class Diagram:
    def __init__(self):
        self._file=""
        self._engine=engine[0]
        self._format="svg"
        self._output=""
        self._cache_str=""

    def _GenerateComment(self):
        with open(self._file+".gv","r") as f:
            flist=[]
            for line in f:
                if line[0:2]=="//":
                    flist.append(line[2:])
        self._cache_str="\n".join(flist)

    def _GenerateImage(self):
        shellstr=self._engine+" -T"+self._format+" "+self._file \
                +".gv -o "+self._output+"."+self._format
        os.system(shellstr)

    def __call__(self,filename,engine,format="svg",output=None):
        '''filename has no .gv'''
        if output is None: output=filename
        if filename is not self._file:
            self._file=filename
            self._GenerateComment()
            self._GenerateImage()
        if not (engine==self._engine and format==self._format and output==self._output):
            self._engine=engine
            self._format=format
            self._output=output
            self._GenerateImage()
        return (self._output+"."+format, self._cache_str)

class MySplitter(QtGui.QSplitter):
    def __init__(self,parent=None):
        super(MySplitter, self).__init__(parent)
        self.walk=Enumerator()
        self.getDiagram=Diagram()
        self._engine=0
        self._diagram=QtGui.QGraphicsView(self)
        self._content = QtGui.QLabel(self)
        self._content.setWordWrap(True)
        self._content.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse);

        scrollText = QtGui.QScrollArea(self)
        scrollText.setWidget(self._content)
        scrollText.setWidgetResizable(True)

        self.addWidget(self._diagram)
        self.addWidget(scrollText)

    def next(self):
        self._updateContent(self.walk(1))
    def prev(self):
        self._updateContent(self.walk(-1))
    def refresh(self):
        self._updateContent(self.walk(0))
    def reset(self, fname=""):
        if fname=="":
            self._diagram.resetTransform()
        else:
            try:
                self._updateContent(self.walk.jumpto(fname))
            except ValueError:
                QtGui.QMessageBox.warning(self,"Error","Diagram file "+fname+".gv does not exist!")
    def zoomin(self):
        self._diagram.scale(1.2,1.2)
    def zoomout(self):
        self._diagram.scale(1.0/1.2,1.0/1.2)
    def retry(self):
        self._engine += 1
        if self._engine>=len(engine):
                self._engine=0
        self.refresh()
    def save(self):
        eng=engine[self._engine]
        output=self.walk(0, IsFullPath=False)
        if not os.path.exists(diagram_export_path):
            os.makedirs(diagram_export_path)
        output=os.path.join(diagram_export_path, output+"_"+eng)
        self.window().statusBar().showMessage("Saving diagram "+output+" with ("+eng+") engine")
        self.getDiagram(self.walk(0),eng,"eps",output)
        with open(output+".txt","w") as f:
            f.write(self._content.text())
        self.window().statusBar().showMessage("diagram "+output+" with ("+eng+") engine saved!")
        
    def _updateContent(self, fname):
        imagefile, comment=self.getDiagram(fname,engine[self._engine])
        svg = QtSvg.QGraphicsSvgItem(imagefile)
        scale=min(self._diagram.width()/svg.boundingRect().width(),
                  self._diagram.height()/svg.boundingRect().height())
        svg.scale(scale,scale)
        scene = QtGui.QGraphicsScene(self)
        scene.addItem(svg)
        self._diagram.setScene(scene)
        self._content.setText(comment)
        self.window().setWindowTitle(fname+".gv with engine "+engine[self._engine]+")")

class Window(QtGui.QMainWindow):
    def __init__(self):
        super(Window, self).__init__()
        self.qsplit=MySplitter(self)
        self.setCentralWidget(self.qsplit)
        self.initToolBar()
        self.initStatusBar()
        self.showMaximized()
        self.qsplit.moveSplitter(self.qsplit.width()/2, 1)
        self.qsplit.refresh()

    def initToolBar(self):      
        prevAction = QtGui.QAction(QtGui.QIcon(icon+"arrow-left.png"), 'Previous', self)
        prevAction.setShortcuts(prev_keys.values())
        prevAction.triggered.connect(self.qsplit.prev)
        nextAction = QtGui.QAction(QtGui.QIcon(icon+"arrow-right.png"), 'Next', self)
        nextAction.setShortcuts(next_keys.values())
        nextAction.triggered.connect(self.qsplit.next)
        retryAction = QtGui.QAction(QtGui.QIcon(icon+"swap.png"), 'Feel lucky', self)
        retryAction.setShortcuts(retry_keys.values())
        retryAction.triggered.connect(self.qsplit.retry)
        saveAction = QtGui.QAction(QtGui.QIcon(icon+"stiffy.png"), 'Save', self)
        saveAction.setShortcuts(save_keys.values())
        saveAction.triggered.connect(self.qsplit.save)
        zoominAction = QtGui.QAction(QtGui.QIcon(icon+"zoom-in-2.png"), 'Zoom In', self)
        zoominAction.setShortcuts(zoomin_keys.values())
        zoominAction.triggered.connect(self.qsplit.zoomin)
        zoomoutAction = QtGui.QAction(QtGui.QIcon(icon+"zoom-out-2.png"), 'Zoom Out', self)
        zoomoutAction.setShortcuts(zoomout_keys.values())
        zoomoutAction.triggered.connect(self.qsplit.zoomout)
        resetAction = QtGui.QAction(QtGui.QIcon(icon+"refresh.png"), 'Reset', self)
        resetAction.setShortcuts(reset_keys.values())
        resetAction.triggered.connect(self.qsplit.reset)
        helpAction = QtGui.QAction(QtGui.QIcon(icon+"speech-bubble-left-2.png"), 'Help', self)
        helpAction.triggered.connect(self.help)

        self.toolbar=QtGui.QToolBar()
        self.toolbar.addSeparator()
        self.toolbar.setIconSize(QtCore.QSize(24,36))
        self.toolbar.addAction(prevAction)
        self.toolbar.addSeparator()
        self.toolbar.addAction(nextAction)
        self.toolbar.addSeparator()
        self.toolbar.addAction(retryAction)
        self.toolbar.addSeparator()
        self.toolbar.addAction(saveAction)
        self.toolbar.addSeparator()
        self.toolbar.addAction(zoominAction)
        self.toolbar.addSeparator()
        self.toolbar.addAction(zoomoutAction)
        self.toolbar.addSeparator()
        self.toolbar.addAction(resetAction)
        self.toolbar.addSeparator()
        self.toolbar.addAction(helpAction)
        self.toolbar.setToolButtonStyle(QtCore.Qt.ToolButtonTextUnderIcon)
        self.addToolBar(QtCore.Qt.RightToolBarArea,self.toolbar)
    
    def initStatusBar(self):
        portal=QtGui.QLineEdit(self)
        portal.setPlaceholderText("Jump to a given diagram file name ...")
        self.connect(portal,QtCore.SIGNAL("returnPressed()"),
                lambda : self.qsplit.reset(portal.text()))
        self.statusBar().addPermanentWidget(portal)

    def help(self):
        QtGui.QMessageBox.information(self,"Shortcuts Help",
                "Prev diagram: "+str(prev_keys.keys())+"\n"+
                "Next diagram: "+str(next_keys.keys())+"\n"+
                "Feel  lucky : "+str(retry_keys.keys())+"\n"+
                "Save diagram: "+str(save_keys.keys())+"\n"+
                "Zoom   in   : "+str(zoomin_keys.keys())+"\n"+
                "Zoom   out  : "+str(zoomout_keys.keys())+"\n"+
                "Reset  view  : "+str(reset_keys.keys()))

if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    ex = Window()
    sys.exit(app.exec_())
