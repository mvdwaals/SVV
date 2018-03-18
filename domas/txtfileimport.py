import os
from dataconst import *

def times():
    results = []
    times = []
    f = open('time.txt',"r")
    for line in f:
        results.append(line.split())

    for i in range(len(results)):
        for j in range(len(results[i])):
            times.append(float(results[i][j]))

    return times

def fuelusedleft():
    results = []
    fuelusedleft = []
    g = open('Fuelusedleft.txt',"r")
    for line in g:
        results.append(line.split())

    for i in range(len(results)):
        fuelusedleft.append(float(results[i][0])*lb_kg)
        
    return fuelusedleft

def fuelusedright():
    results = []
    fuelusedright = []
    h= open('Fuelusedright.txt',"r")
    for line in h:
        results.append(line.split())

    for i in range(len(results)):
        fuelusedright.append(float(results[i][0])*lb_kg)

    return fuelusedright
