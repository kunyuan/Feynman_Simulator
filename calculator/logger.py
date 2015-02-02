#!/usr/bin/env python
import logging
import signal
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 
#add parentdir into PYTHONPATH, where IO module can be found
import IO
workspace = parentdir

log = logging.getLogger()
log.setLevel(logging.INFO)
ch = logging.StreamHandler(sys.stdout)
fh = logging.FileHandler(workspace+'project.log')

ch.setLevel(logging.INFO)
fh.setLevel(logging.INFO)
formatter = logging.Formatter(fmt="[calc][%(asctime)s][%(levelname)s]:\n%(message)s",
        datefmt='%y/%m/%d %H:%M:%S')
ch.setFormatter(formatter)
fh.setFormatter(formatter)
log.addHandler(ch)
log.addHandler(fh)

def Assert(condition, info):
    if not condition:
        log.error(info)
        raise AssertionError

def Abort(info):
    log.error(info)
    raise AssertionError

class DelayedKeyboardInterrupt(object):
    def __enter__(self):
        print "enter sigint"
        self.signal_received = False
        self.old_handler = signal.getsignal(signal.SIGINT)
        signal.signal(signal.SIGINT, self.handler)

    def handler(self, signal, frame):
        self.signal_received = (signal, frame)
        log.info('SIGINT received. Delaying KeyboardInterrupt...')

    def __exit__(self, type, value, traceback):
        print "leave sigint"
        signal.signal(signal.SIGINT, self.old_handler)
        if self.signal_received:
            self.old_handler(*self.signal_received)
