#!/usr/bin/env python 

import numpy as np 
import sys

try : 
    file = sys.argv[1]
    die78 = sys.argv[2] 
    die1 = sys.argv[3]
except : 
    print "Usage %s <pqr file> <die78> <die1>"%sys.argv[0]
    sys.exit() 


dummyIndex = [] 
dummyCoords = [] 
with open(file) as f: 
    index = 0 
    for line in f.readlines() : 
        if not line.split()[0] == "ATOM" : continue 
        if line.split()[2] == 'DUM' : 
            dummyIndex.append(index) 
            dummyCoords.append(line.split()[5:8]) 
        index += 1 

bV = np.array(dummyCoords[1],dtype=float) - np.array(dummyCoords[0],dtype=float) ##bond vector between two dummy atoms

solventPotentials = np.genfromtxt(die78,comments='#')[dummyIndex] 
vacuumPotentials = np.genfromtxt(die1,comments='#')[dummyIndex] 

rxnPotential = solventPotentials - vacuumPotentials 

dV = rxnPotential[1] - rxnPotential[0]
dR = np.sqrt(np.dot(bV, bV) ) 

print -dV/dR     ##Force = - gradient of potential 


