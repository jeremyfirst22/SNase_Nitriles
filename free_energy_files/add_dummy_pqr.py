#!/usr/bin/env python 

import numpy as np
import sys 
import os

try : 
    fileName = sys.argv[1] 
    sep = float(sys.argv[2] ) 
    resName = sys.argv[3]  
    begAtom = sys.argv[4]  
    endAtom = sys.argv[5]  
except : 
    print "Usage: %s < pqr file > < sep (angstroms) of dummy atoms > < residue name > < atom1 of bond vector > < atom2 of bond vector >"%sys.argv[0] 
    exit() 
    
if not os.path.isfile(fileName) : 
    print "ERROR: %s not found!" %fileName 
    exit() 

with open(fileName) as f : 
    fileData = f.readlines() 

lineData = [] 
for line in fileData : 
    lineData.append(line.split() ) 

for index,line in enumerate(lineData) : 
    if line[0] != 'ATOM' : continue 
    if line[3] != resName : continue 
    if line[2] == begAtom : a1Line = index 
    if line[2] == endAtom : a2Line = index 

try : 
    a1Line
    a2Line
except NameError : 
    print "ERROR: Could not find requested atoms" 
    exit() 

##extract coordinates from fields 5,6,7
a1 = np.array(lineData[a1Line][5:8],dtype='f') 
a2 = np.array(lineData[a2Line][5:8],dtype='f') 

## Bond vector (from a1 to a2), and bond length
bV = a2 - a1 
bL = np.sqrt(np.dot(bV, bV)) 
normbV = bV / bL

##Make sure the requested dummy atoms makes sense
assert sep < bL / 2, "Seperation of dummy atoms is too large for the bond vector defined" 

## move from a1 down bond vector to first dummy atom 
##
##  a1 ------ dum1 ----- dum2 ----- a2 
##  |           |_________|         | 
##  |___________|   sep   |_________|
##  |(bL - sep)/2        (bL-sep)/2 |
##  |all in direction of bond Vector| 
##  |_______________________________|
##            bond length (bL)          

dum1 = normbV * (bL - sep)/2 + a1 
dum2 = normbV * (sep )+ dum1
end = normbV * (bL - sep)/2 + dum2

assert end.all() == a2.all(), "ERROR: Bond vector not recovered" 

record1=['ATOM','0','DUM', resName,'0',str(np.round(dum1[0],decimals=3)),str(np.round(dum1[1],decimals=3)),str(np.round(dum1[2],decimals=3)), '0.00', '0.00'] 
record2=['ATOM','0','DUM', resName,'0',str(np.round(dum2[0],decimals=3)),str(np.round(dum2[1],decimals=3)),str(np.round(dum2[2],decimals=3)), '0.00', '0.00'] 

#### Now print out file, adding dummy atom records
def print_remark(record) : 
    print "%-6s"%record[0], 
    print "%3s"%record[1], 
    for item in record[2:] : 
        print item, 
    print 
def print_atom(record) : 
    print "%-6s"%record[0], 
    print "%4s"%record[1],   ## atom number
    if len(record[2])== 4 :        ### This is annoying, but atom names 
        print "%4s"%record[2],     ## get an extra space for 4 character
    else :                         ## names. 
        print " %-3s"%record[2], 
    print "%3s"%record[3],  ## residue name 
    print "%5s"%record[4],   ## residue number
    print "  ", 
    print "%8s%8s%8s"%(record[5],record[6],record[7]), ##coordinates
    print "%7s"%record[8], ##charge
    print "%6s"%record[9] ##VDW radius 
def print_else(record) : 
    for item in record : 
        print item, 
    print 
    
for record in lineData : 
    if record[0] == 'REMARK' : print_remark(record) 
    
for record in record1, record2 : 
    print_atom(record)     

for record in lineData : 
    if record[0] == 'ATOM' : print_atom(record) 

for record in lineData : 
    if record[0] != 'ATOM' and record[0] != 'REMARK' : 
        print_else(record) 



   











