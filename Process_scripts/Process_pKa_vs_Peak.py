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

left, right = 0.1, 0.75
bottom, top = 0.06, 0.98
hspace=0.35

fig, axarr = plt.subplots(3,1,figsize=(3.25,5)) 
fig.subplots_adjust(left=left, bottom=bottom,right=right,top=top,hspace=hspace)
fig.text((right-left)/2+left,top-(top-bottom)/3,           r"Nitrile frequency (cm$^{-1}$)", ha='center', va='bottom')
fig.text((right-left)/2+left,bottom+(top-bottom)/3-0.03 ,           r"Nitrile frequency (cm$^{-1}$)", ha='center', va='bottom')
fig.text((right-left)/2+left,0.01,           r"Lysine p$K_a$",ha='center', va='bottom')
fig.text(0.01,(top-bottom)/2+bottom,         r"Glutatamate p$K_a$",ha='left',va='center',rotation='vertical')
fig.text(0.01,top-(top-bottom-hspace/2)/6          ,         r"Lysine p$K_a$",ha='left',va='center',rotation='vertical')
fig.text(0.01,bottom+(top-bottom-hspace/2)/6          ,         r"Glutamate p$K_a$",ha='left',va='center',rotation='vertical')

glnDict = {}
lysDict = {}
with open('Exp_data/pKas.dat') as f :
    for line in f :
        if not line.startswith('#') :
            key = line.split()[0]
            pKa  = float(line.split()[1] )
            lysDict[key] = pKa
            pKa  = float(line.split()[2] )
            glnDict[key] = pKa

ax1 = axarr[0]
for molec in molecList : 
    ax1.scatter(absMaxDict[molec][0],lysDict[molec],color=colorDict[molec]) 
    ax1.axvline(2162.5,linestyle='--',color='k')
    ax1.axhline(10.4,linestyle='--',color='k')
    ax1.set_ylim([5,11.5]) 
    ax1.text(0.02,0.95,r"\textsf{A}",va='top',ha='left',transform=ax1.transAxes,fontsize=12)

ax2 = axarr[1]
for molec in molecList : 
    ax2.scatter(absMaxDict[molec][0],glnDict[molec],color=colorDict[molec],label=molec) 
    ax2.axvline(2162.5,linestyle='--',color='k')
    ax2.axhline(4.5,linestyle='--',color='k')
    ax2.text(0.02,0.95,r"\textsf{B}",va='top',ha='left',transform=ax2.transAxes,fontsize=12)

ax2 = axarr[2]
for molec in molecList : 
    ax2.scatter(lysDict[molec],glnDict[molec],color=colorDict[molec]) 
    ax2.axvline(10.4,linestyle='--',color='k')
    ax2.axhline(4.5,linestyle='--',color='k')
    ax2.text(0.02,0.95,r"\textsf{C}",va='top',ha='left',transform=ax2.transAxes,fontsize=12)
    ax2.set_xlim(4.4,11)



fig.legend(loc=(0.75,0.380))#,va='center') #,edgecolor='k',framealpha=1) 
    
fig.savefig('figures/pKa_vs_peak.png',format='png',dpi=500) 




