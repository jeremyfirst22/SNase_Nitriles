import numpy as np 
import matplotlib.pyplot as plt 
from scipy.stats import linregress
import os
from matplotlib import rc_file 

ftlsdata = 'Exp_data/ftls_fits.dat'
absdata = 'Exp_data/abs_data2.dat'

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
"DMSO",
"Water", 
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

peakDict, fwhmDict = {}, {}
with open(absdata) as f :
    for line in f.readlines() :
        if line.startswith('#') : continue
        key,peak,peakError,fwhm,fwhmError = line.split()
        if fwhm == "nan" : continue
#        value = float(value) * -1
        fwhmDict[key] = [float(fwhm),float(fwhmError) ]
        peakDict[key] = [float(peak),float(peakError) ]


left, right = 0.18, 0.75
bottom, top = 0.10, 0.96
hspace = 0.4

fig, axarr = plt.subplots(2,1,sharey='all',sharex='none',figsize=(3.42,3.00))
fig.subplots_adjust(left=left, bottom=bottom,right=right,top=top,hspace=hspace)
fig.text((right-left)/2+left,0.50,           r"Hydrogen bond count",ha='center', va='bottom')
fig.text((right-left)/2+left,0.01,           r"FWHM (cm$^{-1}$)   ",ha='center', va='bottom')
fig.text(0.01,(top-bottom)/2+bottom,r"FTLS (cm$^{-1}$ $^{\circ}$C$^{-1}$ )",ha='left',va='center',rotation='vertical')

ax1, ax2 = axarr 

ftlsAccum, ftlsAccum2, fwhmAccum, dataAccum = [],[],[],[]
for index,molec in enumerate(molecList) : 
    ftls, error = ftlsDict[molec]
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
    
    if not molec == "T62X" : 
        ax1.scatter(data,ftls,color=colorDict[molec],zorder=index,label=molec) 
        ax1.errorbar(data,ftls,yerr=error,color=colorDict[molec],zorder=index) 
    else : 
        ax1.scatter(data,ftls,color=colorDict[molec],zorder=index,label="T62X*") 
        ax1.errorbar(data,ftls,yerr=error,color=colorDict[molec],zorder=index)
        ax1.text(data+500,ftls,"*",color=colorDict[molec],va='bottom',ha='left') 

    fwhm,fwhmError = fwhmDict[molec]

    if not molec in ['DMSO', 'Water'] : 
        ax2.scatter(fwhm, ftls, color=colorDict[molec]) 
        ax2.errorbar(fwhm, ftls, yerr=error,color=colorDict[molec]) 
        fwhmAccum.append(fwhm)
        ftlsAccum2.append(ftls)

    ftlsAccum.append(ftls)
    dataAccum.append(data) 


slope, intercept, r_value, p_value, std_error = linregress(dataAccum,ftlsAccum) 
print "r = %0.3f"%r_value
x = np.linspace(min(dataAccum),max(dataAccum),100) 
ax1.plot(x,slope*x+intercept,'k--') 
ax1.text(0.95,0.82,"$r$ = %0.2f"%r_value,transform=ax1.transAxes,va='top',ha='right') 
ax1.text(0.95,0.95,r"\textsf{A}",va='top',ha='right',transform=ax1.transAxes,fontsize=12)

slope, intercept, r_value, p_value, std_error = linregress(fwhmAccum,ftlsAccum2) 
print "r = %0.3f"%r_value
x = np.linspace(np.min(fwhmAccum),np.max(fwhmAccum),100)
ax2.plot(x,slope*x+intercept,'k--') 
ax2.text(0.95,0.82,"$r$ = %0.2f"%r_value,transform=ax2.transAxes,va='top',ha='right') 
ax2.text(0.95,0.95,r"\textsf{B}",va='top',ha='right',transform=ax2.transAxes,fontsize=12)

fig.legend(loc=(right,0.25) ) 
fig.savefig('figures/combined_ftls_figures.png',format='png',dpi=500) 

        


