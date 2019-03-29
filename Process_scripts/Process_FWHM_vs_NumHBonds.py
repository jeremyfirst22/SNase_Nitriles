import numpy as np 
import matplotlib.pyplot as plt 
from scipy.stats import linregress
import os
from matplotlib import rc_file 

absdata = 'Exp_data/abs_data.dat'

rcFile='rc_files/paper.rc'
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

fwhmDict = {}
with open(absdata) as f : 
    for line in f.readlines() : 
        if line.startswith('#') : continue 
        key,peak,fwhm,peakError = line.split()
        if fwhm == "nan" : continue 
#        value = float(value) * -1 
        fwhmDict[key] = float(fwhm)


left, right = 0.18, 0.75
bottom, top = 0.15, 0.95

fig, ax    = plt.subplots(1,1,sharex='all',figsize=(3.25,2))
fig.subplots_adjust(left=left, bottom=bottom,right=right,top=top)
fig.text((right-left)/2+left,0.01,           r"Hydrogen bond count",ha='center', va='bottom')
fig.text(0.01,(top-bottom)/2+bottom,         r"FWHM (cm$^{-1}$)",ha='left',va='center',rotation='vertical')

fwhmAccum, dataAccum = [],[]
for index,molec in enumerate(molecList) : 
    fwhm = fwhmDict[molec]
    if molec == 'Water' : 
        hbondFile = "CNC/hbond/geometry.xvg"
        if not os.path.isfile(hbondFile) : 
            print "%s not found."%hbondFile 
            continue 
        headlines = 0 
        with open(hbondFile) as f : 
            for line in f.readlines() : 
                if line.startswith('#') or line.startswith('@') : 
                    headlines += 1 
                else : 
                    break 
        data = len(np.genfromtxt(hbondFile,skip_header=headlines) ) 
        data *= 5  ##Water simulation only for 20 ns. Extrapolate to 100 ns to compare. 
    elif not molec == "DMSO" : 
        hbondFile = "SNase_%s/hbond/geometry.xvg"%molec
        if not os.path.isfile(hbondFile) : 
            print "%s not found."%hbondFile 
            continue 
        headlines = 0 
        with open(hbondFile) as f : 
            for line in f.readlines() : 
                if line.startswith('#') or line.startswith('@') : 
                    headlines += 1 
                else : 
                    break 
        if not molec == "T62X" : 
            data = len(np.genfromtxt(hbondFile,skip_header=headlines) ) 

        hbondFile = "SNase_%s/hbond/nw_geometry.xvg"%molec
        if not os.path.isfile(hbondFile) : 
            print "%s not found."%hbondFile 
            continue 
        headlines = 0 
        with open(hbondFile) as f : 
            for line in f.readlines() : 
                if line.startswith('#') or line.startswith('@') : 
                    headlines += 1 
                else : 
                    break 
        data += len(np.genfromtxt(hbondFile,skip_header=headlines) ) 
    else : 
        data = 0 
    
    ax.scatter(data,fwhm,color=colorDict[molec],zorder=index,label=molec) 
    #ax.errorbar(data,fwhm,yerr=error,color=colorDict[molec],zorder=index) 

#    if molec == "A58X" : continue 
    fwhmAccum.append(fwhm)
    dataAccum.append(data) 


slope, intercept, r_value, p_value, std_error = linregress(dataAccum,fwhmAccum) 
print "r = %0.3f"%r_value

x = np.linspace(min(dataAccum),max(dataAccum),100) 
ax.plot(x,slope*x+intercept,'k--') 

#ax.set_xscale('log')  
#ax.set_xlim([10**2,10**5]) 
#ax.set_yscale('log')  
ax.legend(loc=(1.00,0.0))
ax.text(0.95,0.95,"r = %0.2f"%r_value,transform=ax.transAxes,va='top',ha='right') 
fig.savefig('figures/fwhm_vs_numhbonds.png',format='png',dpi=500) 

        


