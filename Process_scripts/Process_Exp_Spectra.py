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
"L25X",
"L38X",
"I92X",
"N118X",
"A58X",
"DMSO",
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

figRows,figCols = 11,1

srfDict, pcfDict = {}, {} 
srfStdDict, pcfStdDict = {}, {} 

left, right = 0.08, 0.95 
bottom, top = 0.08, 0.98 
wspace, hspace = 0.05, 0.0 

fig, axarr = plt.subplots(figRows,figCols,sharex='all',sharey='all',figsize=(2.50,5.5) )
fig.subplots_adjust(left=left,bottom=bottom,right=right,top=top,hspace=hspace,wspace=wspace)
fig.text(left+(right-left)/2,0.03, r"Wavenumbers (cm$^{-1}$)" , ha='center', va='center') 
fig.text(0.05,bottom+(top-bottom)/2, r"Normalized absorbance (a.u.)", ha='center', va='center',rotation='vertical') 

equilTime= 00 
forceDir = "force_calc"
f = open('Exp_data/abs_data2.dat','w') 
f.write("#Construct Mean freq. std  FWHM  std\n") 
for index,molec in enumerate(molecList) :
    if molec == "Water" : index -= 1  ##Water on same row as DMSO

    #ax = axarr[index%figRows,index/figRows]
    ax = axarr[index]
    if molec != "Water" : 
        ax.text(0.02,0.88,"%s"%molec,va='top',ha='left',color=colorDict[molec],transform=ax.transAxes)
    else : 
        ax.text(0.98,0.88,"%s"%molec,va='top',ha='right',color=colorDict[molec],transform=ax.transAxes)
    ######
    ###   Read in experimental spectra 
    ######
    num = molec[1:-1]
    if molec == "Water" : num = "Water" 
    if molec == "DMSO" : num = "DMSO" 
    expFiles = glob.glob('Exp_data/SNase_final_steadystate_ir_replicates/%s/*.out'%num) 

    if not expFiles : 
        print "No exp files found for %s" %molec
        continue 

    peaks, fwhms, expData = [], [], [] 
    print "\t",
    for expFile in expFiles : 
        ##For some reason, the last run of DMSO was run with 1/2 the resolution. 
        ##  This is making sure the other two runs have the same x-axis as the last one. 
        ##  Also A58X run 1 
        if (molec == "DMSO" and expFile.split('_')[-1][0] in ['4','3']) or \
            (molec == "A58X" and expFile.split('_')[-1][0] in ['1']) : 
            temp = np.genfromtxt(expFile)
            temp = temp[1::2]
            expData.append(temp) 
        else : 
            expData.append(np.genfromtxt(expFile) )

        with open(expFile) as fIn : 
            for line in fIn.readlines() : 
                if line.startswith("#") : 
                    #if line.startswith("# Mean vibrational frequency:") : 
                    #    peaks.append(line.split()[4]) 
                    if line.startswith("# Frequency of max absorbance:") : 
                        peaks.append(line.split()[5]) 
                        print "%8.3f  "%(float(line.split()[5])), 
                    if line.startswith("# FWHM:") : 
                        fwhms.append( line.split()[2]) 
                else : break 
    print "\n",


    ######
    ###   Calculate average spectrum from all experimental spectra
    ######
    xs = expData[0][:,0]
    avgData = np.zeros_like(xs) 
    for data in expData : 
        assert data[:,0].all() == xs.all() 
        assert len(data[:,1]) == len(xs) 
        avgData += data[:,1]
        ax.plot(data[:,0], data[:,1],'k-',linewidth=0.5) 
    avgData /= len(expData) 

    outFile = "Exp_data/%s.dat"%(molec) 
    np.savetxt(outFile, np.array([xs,avgData]).T, fmt='%8.3f')  

    avgPeak = np.average(np.array(peaks,dtype=float) ) 
    stdPeak = np.std(np.array(peaks    ,dtype=float) ) 

    avgFwhm = np.average(np.array(fwhms,dtype=float) ) 
    stdFwhm = np.std(np.array(fwhms    ,dtype=float) ) 

    ######
    ###   Calculate standard deviation and mean frequency
    ######
    #avg, std = weighted_avg_and_std(expData[:,0],expData[:,1]) 
    print("%-6s %8.3f %8.3f %8.3f %8.3f\n"%(molec, avgPeak, stdPeak, avgFwhm, stdFwhm)), 
    f.write("%-6s %8.3f %8.3f %8.3f %8.3f\n"%(molec, avgPeak, stdPeak, avgFwhm, stdFwhm)) 


    ######
    ###   Plot experimental spectra    
    ######
    #ax.plot(expData[:,0],expData[:,1],color=colorDict[molec],linestyle='-')
    #ax.axvline(absMaxDict[molec][0],linestyle='--',color=colorDict[molec]) 

    ax.plot(xs,avgData,color=colorDict[molec],linestyle='-')
    ax.axvline(avgPeak ,linestyle='--',color=colorDict[molec]) 
    #ax.axvline(xs[np.argmax(avgData)],linestyle='--',color=colorDict[molec]) 

    ######
    ###   Set plot parameters
    ######
    ax.set_xlim([2140,2180]) 
    ax.set_ylim([-0.05,1.05]) 
    ax.set_xticks([2145,2160,2175])
    ax.set_yticks([])

f.close() 

fig.savefig('figures/spectra.png',format='png' )








    

