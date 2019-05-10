import glob 
import numpy as np 
import matplotlib.pyplot as plt 
import os
from os import sys
from matplotlib.colors import LogNorm
import matplotlib.lines as mlines 
from matplotlib import rc_file
import matplotlib.gridspec as gridspec
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

rcFile='rc_files/paper.rc'
rc_file(rcFile) 

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

figRows,figCols = 5,2
simTime =100 ##ns

equilTime = 0 

#field = 'solvent_rxn_field'
srfDict, pcfDict = {}, {} 
srfStdDict, pcfStdDict = {}, {} 


fig = plt.figure(figsize=(6.5,4))

left, right    = 0.08, 0.95
bottom, top    = 0.10, 0.97
hspace, wspace = 0.35, 0.30

plotArray = gridspec.GridSpec(2,2,wspace=wspace,hspace=hspace,left=left,right=right,bottom=bottom,top=top)
srfArray  = gridspec.GridSpecFromSubplotSpec(figRows,figCols,subplot_spec=plotArray[:,0],wspace=0.1,hspace=0.00)
axSRF = plt.subplot(plotArray[0,1])
axPCF = plt.subplot(plotArray[1,1])

fig.text(0.01,top-(top-bottom/2)/2, r"SRF calculated with RF method in APBS $(\frac{k_B T}{e^{-}\AA})$",ha='left',va='center',rotation    ='vertical')
fig.text(left +(right-left-wspace/2)/4 ,0.01, r"Time (ns)",ha='center',va='bottom')
fig.text(right-(right-left-wspace/2)/4 ,0.01, r"PCF -- MD force field $(\frac{k_B T}{e^{-}\AA})$",ha='center',va='bottom')
fig.text(right-(right-left-wspace/2)/4 ,0.51, r"SRF -- MD force field $(\frac{k_B T}{e^{-}\AA})$",ha='center',va='bottom')
fig.text(0.51,   top-(top-bottom-hspace/2)/4, r"SRF -- APBS $(\frac{k_B T}{e^{-}\AA})$",ha='left',va='center',rotation='vertical')
fig.text(0.51,bottom+(top-bottom-hspace/2)/4, r"PCF -- APBS $(\frac{k_B T}{e^{-}\AA})$",ha='left',va='center',rotation='vertical')

fig.text(0.01,0.98,r"\textsf{A}",va='top',ha='left',fontsize=12)
fig.text(0.49,0.98,r"\textsf{B}",va='top',ha='left',fontsize=12)
fig.text(0.49,0.49,r"\textsf{C}",va='top',ha='left',fontsize=12)


gmxSRFs,gmxPCFs = [],[]
apbsSRFs,apbsPCFs = [],[]
for index,molec in enumerate(molecList) :
    axTime = plt.Subplot(fig,srfArray[index%figRows,index/figRows])
    fig.add_subplot(axTime,sharex='all',sharey='all')
    axTime.text(0.96,0.93,"%s"%molec,va='top',ha='right' \
                ,color=colorDict[molec],transform=axTime.transAxes)
    if not index%figRows + 1 == figRows :
           axTime.tick_params(axis="x",labelbottom=False)
    if index/figRows + 1 == figCols :
           axTime.tick_params(axis="y",labelleft=False)
    axTime.set_ylim([10,-20]) 


    #continue 


    timeFile = "SNase_%s/APBS_fixed/rxn_field.out"%(molec) 

    if not os.path.isfile(timeFile) : 
        print "%s not found. Skipping"%timeFile
        continue 
    
    data = np.genfromtxt(timeFile)
    data = data[int(float(equilTime)/simTime *len(data)):]
    time = np.linspace(0,simTime,len(data) ) 

    axTime.scatter(time,data,s=0.5,color=colorDict[molec],marker='.')  

    meanField = np.average(data)
    stdField = np.std(data)

    gmxSRFfile = "SNase_%s/force_calc_ca/SNase_%s.solvent_rxn_field.projected.xvg"%(molec,molec) 
    gmxData = np.genfromtxt(gmxSRFfile)
    gmxSRF = np.average(gmxData)
    gmxSRFstd = np.std(gmxData)

    axSRF.scatter(gmxSRF, meanField, color=colorDict[molec],zorder=index,label=molec)
    axSRF.errorbar(gmxSRF, meanField, xerr=gmxSRFstd, yerr=stdField, color=colorDict[molec],zorder=index) 
    gmxSRFs.append(gmxSRF) 
    apbsSRFs.append(meanField) 

    gmxPCFfile = "SNase_%s/force_calc_ca/SNase_%s.protein_field.projected.xvg"%(molec,molec) 
    gmxData = np.genfromtxt(gmxPCFfile)
    gmxPCF = np.average(gmxData)
    gmxPCFstd = np.std(gmxData)

    apbsPCFfile = "SNase_%s/APBS_fixed/coloumb_field.out"%(molec) 
    apbsData = np.genfromtxt(apbsPCFfile)
    apbsPCF = np.average(apbsData)
    apbsPCFstd = np.std(apbsData) 

    axPCF.scatter(gmxPCF, apbsPCF,color=colorDict[molec],zorder=index,label=molec)
    axPCF.errorbar(gmxPCF, apbsPCF, xerr=gmxPCFstd, yerr=apbsPCFstd, color=colorDict[molec],zorder=index) 
    gmxPCFs.append(gmxPCF) 
    apbsPCFs.append(apbsPCF) 

slope, intercept, r_value, p_value, std_error = linregress(gmxSRFs,apbsSRFs)
x = np.linspace(np.min(gmxSRFs),np.max(gmxSRFs),100) 
axSRF.plot(x,slope*x+intercept,'k--') 
#axSRF.plot(x,x,'b--') 
axSRF.text(0.05,0.95,"$r$ = %.2f"%r_value,ha='left',va='top',transform=axSRF.transAxes) 
print slope, intercept

slope, intercept, r_value, p_value, std_error = linregress(gmxPCFs,apbsPCFs)
x = np.linspace(np.min(gmxPCFs),np.max(gmxPCFs),100) 
axPCF.plot(x,slope*x+intercept,'k--') 
#axPCF.plot(x,x,'b--') 
axPCF.text(0.05,0.95,"$r$ = %.2f"%r_value,ha='left',va='top',transform=axPCF.transAxes) 
print slope, intercept

fig.savefig('figures/combined_APBS_figure.png',format='png') 

plt.close(fig)

