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
#"A90X",
"V66X",
#"A109X",
#"T62X",
#"V23X",
#"L25X",
"L38X",
"I92X" #,
#"N118X",
#"A58X" 
]

absMaxDict = {}  
with open('Exp_data/abs_data2.dat') as f : 
    for line in f : 
        if not line.startswith('#') : 
            key = line.split()[0] 
            absMax, FWHM, error = line.split()[1:4] 
            absMax, FWHM, error = float(absMax), float(FWHM), float(error) 
            absMaxDict[key] = [absMax,FWHM,error]

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

def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return (average, np.sqrt(variance))


#rcFile='rc_files/paper.rc'
#rc_file(rcFile) 
#plt.xkcd() 

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

figRows,figCols = 1,1
simTime =100 ##ns

srfDict, pcfDict = {}, {} 
srfStdDict, pcfStdDict = {}, {} 


left, right = 0.05, 0.95 
bottom, top = 0.05, 0.95 
wspace, hspace = 0.05, 0.0 


fig, axarr = plt.subplots(figRows,figCols,sharex='all',sharey='all',figsize=(1.00,1.00) )
fig.subplots_adjust(left=left,bottom=bottom,right=right,top=top,hspace=hspace,wspace=wspace)
#fig.text(left+(right-left)/2,0.03, r"Wavenumbers (cm$^{-1}$)" , ha='center', va='center') 
#fig.text(0.05,bottom+(top-bottom)/2, r"Normalized absorbance (a.u.)", ha='center', va='center',rotation='vertical') 

equilTime= 00 
forceDir = "force_calc"
#f = open('Exp_data/abs_data2.dat','w') 
#f.write("#Construct Abs Max  STD   \n") 
for index,molec in enumerate(molecList) :
    ax = axarr #[index%figRows,index/figRows]
#    ax.text(0.02,0.88,"%s"%molec,va='top',ha='left',color=colorDict[molec],transform=ax.transAxes)
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
    #avg, std = weighted_avg_and_std(expData[:,0],expData[:,1]) 
    #print "%6s\t%4.2f\t%2.3f"%(molec, avg, std)
    #f.write("%6s\t%4.2f\t%2.3f\t0.0000\n"%(molec, avg, std)) 


    ######
    ###   Plot experimental spectra    
    ######
    ax.plot(expData[:,0],expData[:,1],color=colorDict[molec],linestyle='-')
    #ax.axvline(absMaxDict[molec][0],linestyle='--',color=colorDict[molec]) 

    ######
    ###   Set plot parameters
    ######
    ax.set_xlim([2150,2172]) 
    ax.set_ylim([-0.02,1.05]) 
#    ax.set_xticks([2145,2160,2175])
    ax.set_yticks([])
    ax.set_xticks([])

#f.close() 

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

fig.savefig('figures/TOC_spectra.png',format='png' ,dpi=500)
fig.clf()



fig, axarr = plt.subplots(figRows,figCols,sharex='all',sharey='all',figsize=(1.00,1.00) )
fig.subplots_adjust(left=left,bottom=bottom,right=right,top=top,hspace=hspace,wspace=wspace)
for index,molec in enumerate(molecList) :
    ax = axarr #[index%figRows,index/figRows]
    ax.bar(index,glnDict[molec],color=colorDict[molec],width=0.5) 

    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_xlim([-1.0,3.0])

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

fig.savefig('figures/TOC_pKa.png',format='png' ,dpi=500)
fig.clf()








    

