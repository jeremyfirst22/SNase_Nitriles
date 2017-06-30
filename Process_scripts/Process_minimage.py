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

datafiles = glob.glob('*/minimage/mindist.xvg') 

index=0
for datafile in datafiles : 
    print index, index/figCols, index%figRows
    ax = axarr[index/figCols,index%figCols]

    print datafile
    try : 
        data = np.genfromtxt(datafile,skip_header=27) 
    except IOError : 
        print "No file found for %s %s"%(state,solvent)  
        continue

    data[:,0] = data[:,0] / 1000 
    data[:,1] = data[:,1] * 10 # nm -> Angstroms

    ax.scatter(data[:,0],data[:,1],s=0.1,color='b') 
    ax.set_title(datafile.split('/')[0]) 
    ax.set_xlim(0,50) 
    ax.set_ylim(10,50.0) 

    index +=1
fig.savefig('figures/minimage.png',format='png') 
