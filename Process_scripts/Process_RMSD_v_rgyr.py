import glob 
import numpy as np 
import matplotlib.pyplot as plt 
import os
from os import sys
from matplotlib.colors import LogNorm
import matplotlib.lines as mlines 
import matplotlib.cm as cm 

figCols=2
figRows=3

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0) 
fig.text(0.5,0.04, r"R$_g$ ($\AA$)", ha='center', va='center') 
fig.text(0.08,0.5, r"RMSD ($\AA$)", ha='center', va='center',rotation='vertical') 

try : 
    binSize = int(sys.argv[1] ) 
except IndexError : 
    print "Usage: %s < binSize >"%(sys.argv[0]) 
    print "Setting binSize to default of 100" 
    binSize = 100 



datafiles = glob.glob('*/rgyr/gyrate.xvg') 

index=0
for datafile in datafiles : 
    molec = datafile.split('/')[0]
    datafile2 = "%s/rmsd/without_ter.xvg"%molec
    print index, index/figCols, index%figRows
    ax = axarr[index/figCols,index%figCols]

    print datafile
    try : 
        data = np.genfromtxt(datafile,skip_header=25) 
    except IOError : 
        print "No gryate file found for %s"%(molec) 
        continue
    try : 
        data2= np.genfromtxt(datafile2,skip_header=16) 
    except IOError : 
        print "No RMSD file found for %s"%(molec) 
        continue

    x = data[:,1] * 10 # nm -> Angstroms
    y = data2[:,1] * 10 
    print len(x), len(y) 
    assert len(x) == len(y) 

    while len(x) % binSize != 0 : 
        x = x[:-1]
        y = y[:-1]
    assert len(x) % binSize == 0 
    assert len(y) % binSize == 0 
    
    xs = np.mean(x.reshape(-1,binSize),axis=1) 
    ys = np.mean(y.reshape(-1,binSize),axis=1) 

    colors = cm.brg(np.linspace(0,1,len(ys)) ) 

    for x,y,c in zip(xs,ys,colors) : 
        ax.scatter(x,y,s=0.1,color=c) 

    ax.set_title(molec) 
    ax.set_xlim(15,17) 
    ax.set_ylim(0,2.5) 

    index +=1
fig.savefig('figures/rmsd_v_rgyr.png',format='png') 
