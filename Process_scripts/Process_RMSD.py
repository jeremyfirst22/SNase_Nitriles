import glob 
import numpy as np 
import matplotlib.pyplot as plt 
import os
from os import sys
from matplotlib.colors import LogNorm
import matplotlib.lines as mlines 
from matplotlib import rc_file

figCols=2
figRows=5

rcFile = 'rc_files/paper.rc'
rc_file(rcFile) 

colorDict = {
        "V23X":'chartreuse',
        "L25X":'g',
        "L38X":'c',
        "A58X":'m',
        "T62X":'gold',
        "V66X":'r',
        "A90X":'darkred',
        "I92X":'b',
        "A109X":'darkorange',
        "N118X":'indigo',
        "DMSO": 'k',
        "Water": 'grey'
}

molecList = [
"A90X",
"V66X",
"A109X",
"T62X",
"V23X",
"L25X",
"L38X",
"I92X",
"N118X",
"A58X"
]

if not os.path.isdir('figures') : 
    os.mkdir('figures') 


datafile = 'SNase_WT/rmsd/crystal.xvg'
if not os.path.isfile(datafile) :
    print "No WT found. Quitting. " 
    sys.exit() 
headlines = 0 
with open(datafile) as f: 
    for line in f.readlines() : 
        if line.startswith('#') or line.startswith('@') : 
            headlines += 1 
        else : 
            break 
try : 
    dataWT= np.genfromtxt(datafile,skip_header=headlines) 
except ValueError : 
    print "Error importing %s"%(molec) 
    sys.exit() 
dataWT[:,0] /=  1000  ##ps -> ns
dataWT[:,1] *=  10    ##nm -> Angstrom

left, right = 0.15, 0.99
bottom, top = 0.08, 0.92
wspace,hspace = 0.1,0.3

fig, axarr = plt.subplots(figRows,figCols,sharex='all',sharey='all',figsize=(3.25,5) )
fig.subplots_adjust(wspace=wspace,hspace=hspace,left=left,right=right,top=top,bottom=bottom)
fig.text(left+(right-left)/2,0.01, r"Time (ns)", ha='center', va='bottom')
fig.text(0.01,bottom+(top-bottom)/2, r"RMSD ($\AA$)", ha='left', va='center',rotation='vertical')

maxResidual = 0 
for index,molec in enumerate(molecList) : 
    ax = axarr[index%figRows,index/figRows]

    ax.text(0.50,1.00,"%s"%molec,va='bottom',ha='center',color=colorDict[molec],transform=ax.transAxes)
    #ax.set_xlim(0,100) 

    datafile = 'SNase_%s/rmsd/crystal.xvg'%(molec) 

    if not os.path.isfile(datafile) : 
        print "No file found for %s"%(molec)
        continue 

    headlines = 0 
    with open(datafile) as f: 
        for line in f.readlines() : 
            if line.startswith('#') or line.startswith('@') : 
                headlines += 1 
            else : 
                break 
    try : 
        data1 = np.genfromtxt(datafile,skip_header=headlines) 
    except ValueError : 
        print "Error importing %s"%(molec) 
        continue 

    data1[:,0] /=  1000  ##ps -> ns
    data1[:,1] *= 10  ##nm -> Angstrom

    ax.scatter(dataWT[:,0],dataWT[:,1],color='k',s=0.01) 
    ax.scatter(data1[:,0],data1[:,1],color=colorDict[molec],s=0.01) 

    residual = np.abs(dataWT[2500:,1] - data1[2500:,1] ) 
    if np.max(residual) > maxResidual : maxResidual = np.max(residual) 

    #ax.set_title(mol,color=colorDict[mol]) 
    #ax.set_ylim(0,.35)     


print "Max residual = %5.3f"%maxResidual
fig.savefig('figures/rmsd.png',format='png') 
plt.close() 

