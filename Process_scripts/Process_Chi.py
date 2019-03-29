import glob 
import numpy as np 
import matplotlib.pyplot as plt 
import os
from os import sys
from matplotlib.colors import LogNorm
import matplotlib.lines as mlines 
import matplotlib.gridspec as gridspec
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

left, right = 0.07, 0.90
bottom, top = 0.08, 0.92
wspace = 0.3

fig = plt.figure(figsize=(6.5,5))

fig.text(left +(right-left-wspace/2)/4,0.01, r"Time (ns)                 ",ha='center',va='bottom')
fig.text(right-(right-left-wspace/2)/4,0.01, r"Time (ns)                 ",ha='center',va='bottom')
fig.text(0.01,bottom+(top-bottom)/2         , r"$\chi_1$ (deg)",ha='left',va='center',rotation='vertical')
fig.text(0.55-wspace/4,bottom+(top-bottom)/2, r"$\chi_2$ (deg)",ha='left',va='center',rotation='vertical')

fig.text(0.01,0.99,r"\textsf{A}",va='top',ha='left',fontsize=12)
fig.text(0.47,0.99,r"\textsf{B}",va='top',ha='left',fontsize=12)

plotArray = gridspec.GridSpec(1,2,wspace=wspace,left=left,right=right,bottom=bottom,top=top)
for subplot,dihedral in enumerate(['chi1','chi2']) : 
    axArr = gridspec.GridSpecFromSubplotSpec(figRows,figCols,subplot_spec=plotArray[subplot],wspace=0.1,hspace=0.35)
    for index,molec in enumerate(molecList) : 
        ax = plt.Subplot(fig,axArr[index%figRows,index/figRows])
        fig.add_subplot(ax,sharex='all',sharey='all')

        if not index%figRows + 1 == figRows :
            ax.tick_params(axis="x",labelbottom=False)
        if index/figRows + 1 == figCols :
            ax.tick_params(axis="y",labelleft=False)
        ax.text(0.50,1.00,"%s"%molec,va='bottom',ha='center',color=colorDict[molec],transform=ax.transAxes)

   #     ax.set_xlim(0,100) 
        ax.set_yticks([-180,0,180]) 
        ax.set_ylim([-180,180]) 

   #     continue 
    
        datafile = 'SNase_%s/chi/%s.xvg'%(molec,dihedral) 
    
        if not os.path.isfile(datafile) : 
            print "No file found for %s"%(mol)
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
            print "Error importing %s"%(mol) 
            continue 
    
        data1[:,0] = data1[:,0] / 1000  ##ps -> ns
    
        ax.scatter(data1[:,0],data1[:,1],color=colorDict[molec],s=0.1) 
    
        #ax.set_title(molec,color=colorDict[molec]) 
        #ax.set_ylim(0,.35)     
    
        index +=1
    
fig.savefig('figures/combined_chi.png',format='png') 
plt.close() 

