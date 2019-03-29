#!/usr/bin/env python 

import numpy as np 
import sys 
import os.path

try : 
    fileName = sys.argv[1] 
    resName  = sys.argv[2] 
    bVstart  = sys.argv[3] 
    bVend    = sys.argv[4] 
except : 
    print "Usage: %s < pqr file > < residue name > <atom1 of bond vector > < atom2 of bond vector> "%sys.argv[0] 
    print "\t   Potential will be calculated at the midpoint of the bond vector." 
    exit() 
if not os.path.isfile(fileName) : 
    print "ERROR: %s not found"%fileName 
    exit() 

atoms = [] 
with open(fileName) as f : 
    for line in f.readlines() : 
        items = line.split() 
        if items[0] == 'ATOM' : 
            atoms.append([items[3], items[2], items[8], items[5], items[6], items[7]]) ## Residue, name, Charge, x, y, z 

for index,atom in enumerate(atoms) : 
    if atom[0] == resName :
        if atom[1] == bVstart : 
            atom1 = index 
            atom[2] = 0.00    ##Zero charge on atoms in bond vector (otherwise they dominate the RF) 
        if atom[1] == bVend : 
            atom2 = index 
            atom[2] = 0.00    ##Zero charge on atoms in bond vector (otherwise they dominate the RF) 

try : 
    atom1, atom2
except NameError : 
    print "ERROR: Unable to find %s and %s in residue %s"%(bVstart, bVend, resName) 
    exit() 

##
## Extract charges and coordinates 
##
atoms = np.array(atoms) 
atoms = np.array(atoms[:,2:],dtype='f')  

bV = (atoms[atom2][1:] - atoms[atom1][1:]) ##Bond vector (If arguments are atom 1 = CD and atom 2 = NE, then bV points from CD -> NE
bMP = bV / 2 + atoms[atom1][1:] ##midpoint of bond vector

k = 1 / ( 4 * 3.14 * 8.85 * 10 ** -12) ## k = 1/(4*pi*Eps0)  N m^2 C^-2
k = k * (1.602 * 10 **-19 )**2         ## k : N m^2 / e^2
k = k * 10**10 * 10**10                ## k : N ang^2 / e^2 
k = k / (4.11 * 10**-11 )              ## k : kbT / A / Ang^2 / e^2


##
## Vector coloumbs: 
##                                
##   F = k * q1*q2/ magR^2  *unitR
##   F = k * q1*q2 / magR**2* (vecR / magR)
##   F = k * q1*q1 / magR**3 * vecR
##         
F = 0 
for atom in atoms : 
    q = atom[0] 
    rVec = bMP - atom[1:] 
    r = np.sqrt(np.dot(rVec,rVec))

    F += q * rVec / r**3 
F *= k

#print F   ## Potential for charge e, is identical in units of kbT/eA 
print np.dot(F, bV) 
#print np.sqrt(F[0]**2 + F[1]**2 + F[2]**2 ) 

