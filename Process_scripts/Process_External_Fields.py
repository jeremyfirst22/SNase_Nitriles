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

absMaxDict = {}
with open('Exp_data/abs_data2.dat') as f :
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

figRows,figCols = 5,2
simTime =100 ##ns

equilTime = 0 

#field = 'solvent_rxn_field'
srfDict, pcfDict = {}, {} 
srfStdDict, pcfStdDict = {}, {} 


fig = plt.figure(figsize=(3.25,4))

left, right    = 0.15, 0.75
bottom, top    = 0.10, 0.97
hspace, wspace = 0.35, 0.30

plotArray = gridspec.GridSpec(2,1,wspace=wspace,hspace=hspace,left=left,right=right,bottom=bottom,top=top)
gmxAx = plt.subplot(plotArray[0])
apbsAx = plt.subplot(plotArray[1])

fig.text(0.02,top-(top-bottom-hspace/2)/4, r"$\langle F \rangle$ -- MD force field $(\frac{k_B T}{e^{-}\AA})$",ha='left',va='center',rotation    ='vertical')
fig.text(0.02,bottom+(top-bottom-hspace/2)/4, r"$\langle F \rangle$ -- APBS $(\frac{k_B T}{e^{-}\AA})$",ha='left',va='center',rotation    ='vertical')
fig.text(left +(right-left)/2 ,0.01, r"$\tilde{\nu}$ (cm$^{-1}$)",ha='center',va='bottom')
fig.text(left +(right-left)/2 ,0.51, r"$\tilde{\nu}$ (cm$^{-1}$)",ha='center',va='bottom')

gmxAx.text(0.05,0.95,r"\textsf{A}",va='top',ha='left',fontsize=12,transform=gmxAx.transAxes)
apbsAx.text(0.05,0.95,r"\textsf{B}",va='top',ha='left',fontsize=12,transform=apbsAx.transAxes)

gmxEFs, apbsEFs, absMaxs = [],[],[]
gmxSRFstds, gmxPCFstds = [], [] 
fwhms = [] 
for index,molec in enumerate(molecList) :

    gmxSRFfile = "SNase_%s/force_calc_ca/SNase_%s.solvent_rxn_field.projected.xvg"%(molec,molec) 
    gmxData = np.genfromtxt(gmxSRFfile)
    gmxSRF = np.average(gmxData)
    gmxStdSRF = np.std(gmxData)

    apbsSRFfile = "SNase_%s/APBS_fixed/rxn_field.out"%(molec) 
    apbsData = np.genfromtxt(apbsSRFfile)
    apbsSRF = np.average(apbsData)


    gmxPCFfile = "SNase_%s/force_calc_ca/SNase_%s.protein_field.projected.xvg"%(molec,molec) 
    gmxData = np.genfromtxt(gmxPCFfile)
    gmxPCF = np.average(gmxData)
    gmxStdPCF = np.std(gmxData)

    apbsPCFfile = "SNase_%s/APBS_fixed/coloumb_field.out"%(molec) 
    apbsData = np.genfromtxt(apbsPCFfile)
    apbsPCF = np.average(apbsData)

    gmxEF = gmxSRF + gmxPCF
    apbsEF = apbsSRF + apbsPCF

    gmxAx.scatter(absMaxDict[molec][0],gmxEF,color=colorDict[molec],label=molec) 
    apbsAx.scatter(absMaxDict[molec][0],apbsEF,color=colorDict[molec]) 

    gmxEFs.append(gmxEF)
    apbsEFs.append(apbsEF)
    absMaxs.append(absMaxDict[molec][0]) 

    gmxSRFstds.append(gmxStdSRF)
    gmxPCFstds.append(gmxStdPCF)
    fwhms.append(absMaxDict[molec][1]) 


x = np.linspace(np.min(absMaxs),np.max(absMaxs),100) 

slope, intercept, r_value, p_value, std_error = linregress(absMaxs,gmxEFs)
gmxAx.plot(x,slope*x+intercept,'k--') 
gmxAx.text(0.05,0.05,"r = %.2f"%r_value,ha='left',va='bottom',transform=gmxAx.transAxes) 
print slope, intercept

slope, intercept, r_value, p_value, std_error = linregress(absMaxs,apbsEFs)
apbsAx.plot(x,slope*x+intercept,'k--') 
apbsAx.text(0.05,0.05,"r = %.2f"%r_value,ha='left',va='bottom',transform=apbsAx.transAxes) 
print slope, intercept

fig.legend(loc=(0.75,0.35))
fig.savefig('figures/combined_external_fields_figure.png',format='png') 

plt.close(fig)

fig = plt.figure(figsize=(3.25,4))

plotArray = gridspec.GridSpec(2,1,wspace=wspace,hspace=hspace,left=left,right=right,bottom=bottom,top=top)
srfAx = plt.subplot(plotArray[0])
pcfAx = plt.subplot(plotArray[1])

fig.text(0.02,top-(top-bottom-hspace/2)/4, r"$\sigma_{SRF}$ $(\frac{k_B T}{e^{-}\AA})$",ha='left',va='center',rotation    ='vertical')
fig.text(0.02,bottom+(top-bottom-hspace/2)/4, r"$\sigma_{PCF}$ $(\frac{k_B T}{e^{-}\AA})$",ha='left',va='center',rotation    ='vertical')
fig.text(left +(right-left)/2 ,0.01, r"FWHM (cm$^{-1}$)",ha='center',va='bottom')
fig.text(left +(right-left)/2 ,0.51, r"FWHM (cm$^{-1}$)",ha='center',va='bottom')

srfAx.text(0.05,0.95,r"\textsf{A}",va='top',ha='left',fontsize=12,transform=gmxAx.transAxes)
pcfAx.text(0.05,0.95,r"\textsf{B}",va='top',ha='left',fontsize=12,transform=apbsAx.transAxes)

for index,molec in enumerate(molecList) :
    srfAx.scatter(absMaxDict[molec][1],gmxSRFstds[index],color=colorDict[molec],label=molec) 
    pcfAx.scatter(absMaxDict[molec][1],gmxPCFstds[index],color=colorDict[molec]) 

x = np.linspace(np.min(fwhms),np.max(fwhms),100) 

slope, intercept, r_value, p_value, std_error = linregress(fwhms,gmxSRFstds)
srfAx.plot(x,slope*x+intercept,'k--') 
srfAx.text(0.95,0.05,"r = %.2f"%r_value,ha='right',va='bottom',transform=gmxAx.transAxes) 
print slope, intercept

slope, intercept, r_value, p_value, std_error = linregress(fwhms,gmxPCFstds)
pcfAx.plot(x,slope*x+intercept,'k--') 
pcfAx.text(0.95,0.05,"r = %.2f"%r_value,ha='right',va='bottom',transform=apbsAx.transAxes) 
print slope, intercept

fig.legend(loc=(0.75,0.35))
fig.savefig('figures/combined_standard_deviations.png',format='png') 

