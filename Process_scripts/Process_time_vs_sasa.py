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


for sasa in ['sidechain_area'] : #['cnc_area','sidechain_area', 'thiocyanate','nitrile_area'] : 

    left, right = 0.12, 0.95
    bottom, top = 0.08, 0.92
    wspace,hspace = 0.15,0.3

    fig, axarr = plt.subplots(figRows,figCols,sharex='all',sharey='all',figsize=(3.25,5))
    fig.subplots_adjust(left=left,right=right,top=top,bottom=bottom,hspace=hspace,wspace=wspace) 
    fig.text(left+(right-left)/2,0.01, r"Time (ns)", ha='center', va='bottom')
    fig.text(0.01,bottom+(top-bottom)/2, r"SASA ($\AA^2$)", ha='left', va='center',rotation='vertical')

    for index, molec in enumerate(molecList) : 
        ax = axarr[index%figRows,index/figRows]
        ax.text(0.50,1.00,"%s"%molec,va='bottom',ha='center',color=colorDict[molec],transform=ax.transAxes)
        ax.set_xlim(0,100) 
        ax.set_ylim(0,35.0)

        datafile = 'SNase_%s/sasa/%s.xvg'%(molec,sasa) 

        if not os.path.isfile(datafile) : 
            print "No file found for %s %s"%(molec, sasa) 
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
            print "Error importing %s %s"%(molec, sasa)
            continue 

        data1[:,0] = data1[:,0] / 1000  ##ps -> ns
        data1[:,2] *= 100               ##nm^2 -> Angstrom^2
#       equilTime = len(data1) / 5 
#       data1 = data1[equilTime:] ##Discard first 10ns as equilibration time
        ax.scatter(data1[:,0],data1[:,2],color=colorDict[molec],s=0.1) 
    
    
    
    fig.savefig('figures/sasa_%s_v_time.png'%(sasa.split('_')[0]),format='png') 
    plt.close() 

