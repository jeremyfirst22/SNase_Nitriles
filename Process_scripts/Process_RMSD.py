import glob 
import numpy as np 
import matplotlib.pyplot as plt 
import os
from os import sys
from matplotlib.colors import LogNorm
import matplotlib.lines as mlines 

figCols=2
figRows=3

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0) 
fig.text(0.5,0.04, "Time (ns)", ha='center', va='center') 
fig.text(0.08,0.5, r"Distance to nearest image ($\AA$)", ha='center', va='center',rotation='vertical') 

datafiles = glob.glob('*/rmsd/backbone.xvg') 

index=0
for datafile in datafiles : 
    print index, index/figCols, index%figRows
    ax = axarr[index/figCols,index%figCols]

    molec = datafile.split('/')[0] 
    datafile2 = molec+"/rmsd/without_ter.xvg"
    print datafile2
    try : 
        data2 = np.genfromtxt(datafile2,skip_header=27) 
    except IOError : 
        print "%s not found!"%datafile2

    print datafile
    try : 
        data = np.genfromtxt(datafile,skip_header=27) 
    except IOError : 
        print "No file found for %s %s"%(state,solvent)  
        continue

    data[:,0] = data[:,0] / 1000 
    data[:,1] = data[:,1] * 10 # nm -> Angstroms

    data2[:,0] = data2[:,0] / 1000 
    data2[:,1] = data2[:,1] * 10 # nm -> Angstroms

    ax.scatter(data[:,0],data[:,1],s=0.1,color='b') 
    ax.scatter(data2[:,0],data2[:,1],s=0.1,color='g') 

    ax.set_title(molec) 
    ax.set_xlim(0,50) 
    ax.set_ylim(0,6.0) 

    index +=1
fig.savefig('figures/rmsd.png',format='png') 