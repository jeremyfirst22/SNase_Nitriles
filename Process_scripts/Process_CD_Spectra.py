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

molecList = [
    "V23X", 
    "L25X", 
    "L38X", 
    "T62X"
]



rcFile='rc_files/paper.rc'
rc_file(rcFile) 

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

figRows,figCols = 2,2
simTime =100 ##ns


left, right = 0.17, 0.99 
bottom, top = 0.11, 0.99 
wspace, hspace = 0.05, 0.05

fig, axarr = plt.subplots(figRows,figCols,sharex='all',sharey='all',figsize=(3.25,3) )
fig.subplots_adjust(left=left,bottom=bottom,right=right,top=top,hspace=hspace,wspace=wspace)
fig.text(left+(right-left)/2,0.01, r"Wavenumbers (cm$^{-1}$)" , ha='center', va='bottom') 
fig.text(0.01,bottom+(top-bottom)/2, r"Ellipticity (x10$^3$ mdeg cm$^2$ dmol$^{-1}$)", ha='left', va='center',rotation='vertical') 

expFile = 'Exp_data/cd_spectra.dat'
if not os.path.isfile(expFile) : 
    print "Experimental spectra not found for %s"%molec
    #continue 

expData = np.genfromtxt(expFile,skip_header=1) 
print expData

for index,molec in enumerate(molecList) :
    ax = axarr[index%figRows,index/figRows]
    ax.text(0.95,0.95,"%s"%molec,va='top',ha='right',color=colorDict[molec],transform=ax.transAxes)
    ######
    ###   Read in experimental spectra 
    ######
    x = expData[:,0] 
    y1= expData[:,index+1]
    y2= expData[:,index+1+len(molecList)]

    y1 /= 1000
    y2 /= 1000

    ######
    ###   Calculate standard deviation and mean frequency
    ######
    #avg, std = weighted_avg_and_std(expData[:,0],expData[:,1]) 
    #print "%6s\t%4.2f\t%2.3f"%(molec, avg, std)
    #f.write("%6s\t%4.2f\t%2.3f\t0.0000\n"%(molec, avg, std)) 


    ######
    ###   Plot experimental spectra    
    ######
    ax.plot(x,y1,color=colorDict[molec],linestyle='-')
    ax.plot(x,y2,color=colorDict[molec],linestyle='--')

    ######
    ###   Set plot parameters
    ######
#    ax.set_xlim([2140,2180]) 
#    ax.set_ylim([-0.05,1.05]) 
#    ax.set_xticks([2145,2160,2175])
#    ax.set_yticks([])


fig.savefig('figures/cd_spectra.png',format='png' )








    

