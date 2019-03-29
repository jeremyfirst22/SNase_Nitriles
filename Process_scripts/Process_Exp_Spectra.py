import glob 
import numpy as np 
import matplotlib.pyplot as plt 
import os
from os import sys
from matplotlib.colors import LogNorm
import matplotlib.lines as mlines 
from matplotlib import rc_file
from scipy.stats import linregress

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
        "DMSO":'k',
        "Water":'grey' 
        }
shapeDict = {
        "V23X":'o',
        "L25X":'o',
        "L38X":'s',
        "A58X":'s',
        "T62X":'D',
        "V66X":'D',
        "A90X":'o',
        "I92X":'D',
        "A109X":'P',
        "N118X":'D'
}
molecList = [
"A90X",
"V66X",
"A109X",
"T62X",
"V23X",
"DMSO",
"L25X",
"L38X",
"I92X",
"N118X",
"A58X",
"Water"
]

absMaxDict = {}  
with open('Exp_data/abs_data2.dat') as f : 
    for line in f : 
        if not line.startswith('#') : 
            key = line.split()[0] 
            absMax, FWHM, error = line.split()[1:4] 
            absMax, FWHM, error = float(absMax), float(FWHM), float(error) 
            absMaxDict[key] = [absMax,FWHM,error]


def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return (average, np.sqrt(variance))


rcFile='rc_files/paper.rc'
rc_file(rcFile) 

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

figRows,figCols = 6,2
simTime =100 ##ns

srfDict, pcfDict = {}, {} 
srfStdDict, pcfStdDict = {}, {} 


left, right = 0.08, 0.99 
bottom, top = 0.11, 0.99 
wspace, hspace = 0.05, 0.0 

fig, axarr = plt.subplots(figRows,figCols,sharex='all',sharey='all',figsize=(3.25,3) )
fig.subplots_adjust(left=left,bottom=bottom,right=right,top=top,hspace=hspace,wspace=wspace)
fig.text(left+(right-left)/2,0.03, r"Wavenumbers (cm$^{-1}$)" , ha='center', va='center') 
fig.text(0.05,bottom+(top-bottom)/2, r"Normalized absorbance (a.u.)", ha='center', va='center',rotation='vertical') 

equilTime= 00 
forceDir = "force_calc"
f = open('Exp_data/abs_data2.dat','w') 
f.write("#Construct Abs Max  STD   \n") 
for index,molec in enumerate(molecList) :
    ax = axarr[index%figRows,index/figRows]
    ax.text(0.02,0.88,"%s"%molec,va='top',ha='left',color=colorDict[molec],transform=ax.transAxes)
    ######
    ###   Read in experimental spectra 
    ######
    expFile = 'Exp_data/%s.dat'%molec
    if not os.path.isfile(expFile) : 
        print "Experimental spectra not found for %s"%molec
        continue 
    expData = np.genfromtxt(expFile) 

    ######
    ###   Calculate standard deviation and mean frequency
    ######
    avg, std = weighted_avg_and_std(expData[:,0],expData[:,1]) 
    print "%6s\t%4.2f\t%2.3f"%(molec, avg, std)
    f.write("%6s\t%4.2f\t%2.3f\t0.0000\n"%(molec, avg, std)) 


    ######
    ###   Plot experimental spectra    
    ######
    ax.plot(expData[:,0],expData[:,1],color=colorDict[molec],linestyle='-')
    ax.axvline(absMaxDict[molec][0],linestyle='--',color=colorDict[molec]) 

    ######
    ###   Set plot parameters
    ######
    ax.set_xlim([2140,2180]) 
    ax.set_ylim([-0.05,1.05]) 
    ax.set_xticks([2145,2160,2175])
    ax.set_yticks([])

f.close() 

fig.savefig('figures/spectra.png',format='png' )








    

