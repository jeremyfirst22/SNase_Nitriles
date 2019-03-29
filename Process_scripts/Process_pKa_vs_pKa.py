import glob 
import numpy as np 
import matplotlib.pyplot as plt 
import os
from os import sys
from matplotlib.colors import LogNorm
import matplotlib.lines as mlines 
from matplotlib import rc_file
from scipy.stats import linregress

figCols=1
figRows=1

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
        "N118X":'indigo'
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

rc_file('rc_files/paper.rc') 

absMaxDict = {}  
with open('Exp_data/abs_data.dat') as f : 
    for line in f : 
        if not line.startswith('#') : 
            key = line.split()[0] 
            absMax, FWHM, error = line.split()[1:4] 
            absMax, FWHM, error = float(absMax), float(FWHM), float(error) 
            absMaxDict[key] = [absMax,FWHM,error]


rcFile='rc_files/paper.rc'
rc_file(rcFile) 

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

left, right = 0.1, 0.8
bottom, top = 0.2, 0.98
hspace=0.1

fig, ax    = plt.subplots(1,1,sharex='all',figsize=(3.25,2)) 
fig.subplots_adjust(left=left, bottom=bottom,right=right,top=top,hspace=hspace)
fig.text((right-left)/2+left,0.01,                  r"Lysine p$K_a$", ha='center', va='bottom')
fig.text(0.01,(top-bottom-hspace)/2+bottom,         r"Glutatamate p$K_a$                ",ha='left',va='center',rotation='vertical')

lysDict = {}
gluDict = {}
with open('Exp_data/pKas.dat') as f :
    for line in f :
        if not line.startswith('#') :
            key = line.split()[0]
            pKa  = float(line.split()[1] )
            lysDict[key] = pKa
            pKa  = float(line.split()[2] )
            gluDict[key] = pKa

#ax = axarr[0]
for molec in molecList : 
    ax.scatter(lysDict[molec],gluDict[molec],color=colorDict[molec],label=molec) 
    #ax.axvline(2162.5,linestyle='--',color='k')
    ax.axvline(10.4,linestyle='--',color='k')
    ax.axhline(4.5,linestyle='--',color='k')

fig.legend(loc=(0.78,0.150)) #,edgecolor='k',framealpha=1) 
    
fig.savefig('figures/pKa_vs_pKa.png',format='png',dpi=500) 




