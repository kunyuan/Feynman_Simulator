#!/usr/bin/env python
 
f = open("most_common_object.txt")
 
alltimings = {}
things = {}
 
start = None
for line in f.readlines():
    if line.startswith("#START"):
        alltimings[start] = things
 
        _, _, thetime = line.strip().partition(" ")
        start = float(thetime)
        things = {}
 
    else:
        key, _, value = line.partition(" ")
        things[key] = int(value.strip())
 
 
times = {}
values = {}
 
import datetime
 
for time, things in sorted(alltimings.items()):
    for t, howmany in things.items():
        times.setdefault(t, []).append(datetime.datetime.fromtimestamp(time))
        values.setdefault(t, []).append(howmany)
 
 
import matplotlib.pyplot as plt
 
thing = plt.figure(figsize = (15, 200))
thing.subplots_adjust(hspace = .2)
 
for i, key in enumerate(times):
    ct = times[key]
    cv = values[key]
 
    ax = thing.add_subplot(len(times.keys()), 1, i+1)
    ax.set_title("%s" % key)
    ax.set_ylabel("Number of objects")
 
    ax.plot(ct, cv)
 
plt.savefig("mem.pdf", figsize = (200, 200), ymargin = 100)
