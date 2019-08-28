import numpy as np 
import matplotlib.pyplot as plt 
import math
from scipy.stats import linregress 
from matplotlib import rc_file
import glob as glob 
from sys import exit

#sys.path.insert(0, "/Users/jfirst/
from adjustText import adjust_text

file = 'Exp_data/FTLS.dat'
outfile = 'Exp_data/ftls_fits.dat'

tempList = [] 
ftlsDict = {} 
with open(file) as f : 
    for item in f.readline().split() : 
        if "std" in item : continue 
        if "#" in item : continue 
        if "deg" in item : continue 
        else : 
            tempList.append(float(item)) 

    for line in f.readlines() : 
        molec = line.split()[0]
        peaks = line.split()[1::2]
        stds  = line.split()[2::2]
        peaks = np.array(peaks, dtype=float) 
        stds = np.array(stds, dtype=float) 
        ftlsDict[molec] = [peaks,stds]

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
    #"Water",
    #"DMSO"
]

rc_file('rc_files/paper.rc') 

left, right = 0.15, 0.99
bottom, top = 0.11, 0.99
wspace, hspace = 0.05, 0.0

figRows, figCols = 5,2
fig, axarr = plt.subplots(figRows,figCols,figsize=(3.42,3.00),sharex='none',sharey='all') 
fig.subplots_adjust(left=left, bottom=bottom,right=right,top=top,hspace=hspace, wspace=wspace) 

fig.text((right-left)/2+left,0.01,           r"Temperature ($^{\circ}$C)",ha='center', va='bottom')
fig.text(0.01,(top-bottom)/2+bottom,         r"Nitrile frequency shift (cm$^{-1}$)                ",ha='left',va='center',rotation='vertical')

f = open(outfile,'w') 

texts = [] 
#for i in range(1,len(data[0])) : 
for index, molec in enumerate(molecList) : 
    ax = axarr[index%figRows,index/figRows]
    #ax.set_title(molec) 
    temps = np.array(tempList,dtype=float) 
    num = molec[1:-1]

    p1Accum, p2Accum = [], [] 
    freq5degree1, freq5degree2 = None, None
    for temp in temps : 
        peak1, peak2,maxPeak = [], [], [] 
        datafiles = glob.glob('Exp_data/SNase_final_ftls_replicates/%s/*_%ideg.out'%(num,temp)) 
        
        for file in datafiles : 
            p1, p2 = None, None
            with open(file) as f : 
                for line in f.readlines() : 
                    if line.startswith('#Gaussian #1') : 
                        p1 = line.split()[7]
                    if line.startswith('#Gaussian #2') : 
                        p2 = line.split()[7]
                    if line.startswith('# Frequency of max absorbance:') : 
                        maxPeak.append(line.split()[5]) 
        #    print p1, p2

            if p1 and p2 : 
                if p2 < p1 : 
                    pTemp = p1
                    p1 = p2 
                    p2 = pTemp

            peak1.append(p1) 
            if p2 : 
                peak2.append(p2) 

        peak1 = np.array(peak1,dtype=float) 
        peak2 = np.array(peak2,dtype=float) 

        peak1Avg = np.average(peak1) 
        peak2Avg = np.average(peak2) 

        #if temp == 5 : 
        #    freq5degree1 = peak1Avg
        #    freq5degree2 = peak2Avg
        #if freq5degree1 and freq5degree2 : 
        #    peak1Avg -= freq5degree1
        #    peak2Avg -= freq5degree2
        #else : 
        #    print "Uh oh. 5degrees undefined"
        #    exit() 

        #print np.average(peak1), np.average(peak2)

        #print temp, np.average(peak1) 
        #maxPeak = np.array(maxPeak,dtype=float) 

        #peak1 -= np.average(maxPeak) 
        #peak2 -= np.average(maxPeak) 

        #for peak in peak1 : 
        ax.scatter(temp, peak1Avg,color=colorDict[molec],marker='d') 
        #for peak in peak2 : 
        ax.scatter(temp, peak2Avg,color=colorDict[molec]) 
        p1Accum.append(peak1Avg) 
        p2Accum.append(peak2Avg) 

    ax.set_xticks([5,15,25,35]) 
    #ax.set_yticks([0,-4.0,-8.0]) 
    ax.text(0.95,0.9,molec,color=colorDict[molec],transform=ax.transAxes,ha='right',va='top') 

    p1Accum = np.array(p1Accum, dtype=float) 
    p2Accum = np.array(p2Accum, dtype=float) 

    print molec
    slope, intercept, r_value, p_values, std_error = linregress(temps, p1Accum) 
    print "\t%5.3f\t%5.3f"%(slope,std_error) 
    xs = np.linspace(np.min(temps), np.max(temps), 100) 
    ax.plot(xs, slope*xs+intercept, color=colorDict[molec],linestyle='-') 

    slope, intercept, r_value, p_values, std_error = linregress(temps, p2Accum) 
    print "\t%5.3f\t%5.3f"%(slope,std_error) 
    xs = np.linspace(np.min(temps), np.max(temps), 100) 
    ax.plot(xs, slope*xs+intercept, color=colorDict[molec],linestyle='--') 


fig.savefig('figures/ftls_both_gaussians.png',format='png',dpi=500) 
