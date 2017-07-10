import glob 
import numpy as np 
import matplotlib.pyplot as plt 
import os
from os import sys
from matplotlib.colors import LogNorm
import matplotlib.lines as mlines 
#from matplotlib import rc_file 

figCols=2
figRows=3

#rcFile = 'rc_files/paper.rc' 

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0.1) 
fig.subplots_adjust(hspace=0.25) 
fig.text(0.5,0.04, "Time (ns)", ha='center', va='center') 
fig.text(0.05,0.5, r"Distance to nearest image ($\rm{\AA}$)", ha='center', va='center',rotation='vertical') 

#rc_file(rcFile) 

datafiles = glob.glob('*/minimage/mindist.xvg') 

index=0
for datafile in datafiles : 
    molec = datafile.split('/')[0]
    title = molec.split('_')[1]
    if title[-1] == 'X' : 
        title = title[:-1]
        title = title+'C$_{\\rm{SCN}}$'
    print index, index/figCols, index%figRows
    ax = axarr[index/figCols,index%figCols]

    print datafile
    try : 
        data = np.genfromtxt(datafile,skip_header=27) 
    except IOError : 
        print "No file found for %s %s"%(state,solvent)  
        continue
    except ValueError : 
        print "Trying without last line"
        data = np.genfromtxt(datafile,skip_header=27,skip_footer=1) 

    data[:,0] = data[:,0] / 1000 
    data[:,1] = data[:,1] * 10 # nm -> Angstroms

    ax.scatter(data[:,0],data[:,1],s=0.1,color='b') 
    ax.set_title(title) 
    ax.set_xlim(0,50) 
    ax.set_ylim(10,50.0) 

    #if not index % figCols == 0 : 
    #    xticks = ax.xaxis.get_major_ticks() 
    #    xticks[0].label1.set_visible(False) 

    index +=1
fig.savefig('figures/minimage.png',format='png') 
