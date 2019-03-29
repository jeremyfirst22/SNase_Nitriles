import numpy as np 
import matplotlib.pyplot as plt 
from scipy.stats import linregress
import os
from matplotlib import rc_file 

ftlsdata = 'Exp_data/ftls.dat'
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

ftlsDict = {}
with open(ftlsdata) as f : 
    for line in f.readlines() : 
        key,value, error = line.split()
        if value == "nan" : continue 
#        value = float(value) * -1 
        ftlsDict[key] = [float(value), float(error)] 

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
fig.text((right-left)/2+left,0.01,           r"FWHM (cm$^{-1}$)",ha='center', va='bottom')
fig.text(0.01,(top-bottom)/2+bottom,         r"FTLS (cm$^{-1}$ $^{\circ}$C$^{-1}$ )",ha='left',va='center',rotation='vertical')

ftlsAccum, fwhmAccum = [],[]
for index,molec in enumerate(molecList) : 
    ftls, error = ftlsDict[molec]
    fwhm = fwhmDict[molec]
    
    ax.scatter(fwhm,ftls,color=colorDict[molec],zorder=index,label=molec) 
    ax.errorbar(fwhm,ftls,yerr=error,color=colorDict[molec],zorder=index) 

#    if molec == "A58X" : continue 
    ftlsAccum.append(ftls)
    fwhmAccum.append(fwhm) 


slope, intercept, r_value, p_value, std_error = linregress(fwhmAccum,ftlsAccum) 
print "r = %0.3f"%r_value

x = np.linspace(min(fwhmAccum),max(fwhmAccum),100) 
ax.plot(x,slope*x+intercept,'k--') 

#ax.set_xscale('log')  
#ax.set_xlim([10**2,10**5]) 
#ax.set_yscale('log')  
ax.legend(loc=(1.00,0.0))
ax.text(0.95,0.95,"r = %0.2f"%r_value,transform=ax.transAxes,va='top',ha='right') 
fig.savefig('figures/ftls_vs_fwhm.png',format='png',dpi=500) 

        


