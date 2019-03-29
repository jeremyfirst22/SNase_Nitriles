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
        "V23X":'#007F00',
        "L25X":'#4DBEEE',
        "L38X":'#142B8C',
        "A58X":'#BF00BF',
        "T62X":'#77AB30',
        "V66X":'#D95219',
        "A90X":'#A2142F',
        "I92X":'#6666FF',
        "A109X":'#ECB120',
        "V104X":'k',
        "N118X":'#7E2F8E'
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
"A58X"
]

absMaxDict = {}  
with open('Exp_data/abs_data.dat') as f : 
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

srfDict, pcfDict = {}, {} 
srfStdDict, pcfStdDict = {}, {} 


fig, axarr = plt.subplots(figRows,figCols,sharex='all',sharey='all',figsize=(4.3,6) )
fig.subplots_adjust(left=0.10,bottom=0.10,right=0.95,top=0.90,hspace=0.6,wspace=0.2)
fig.text(0.5,0.04, r"Shift due to field (cm$^{-1}$)" , ha='center', va='center') 
fig.text(0.5,0.97, r"Experimental frequency (cm$^{-1}$)" , ha='center', va='center') 
fig.text(0.03,0.5, r"Occurances", ha='center', va='center',rotation='vertical') 

#fig2, ax2 = plt.subplots(1,1,figsize=(4.3,3) )
#fig2.subplots_adjust(left=0.20,bottom=0.13,right=0.95,top=0.95)
#fig2.text(0.5,0.02, r"Shift due to %s (cm$^{-1}$)"%field.replace('_',' '),ha='center', va='center') 
#fig2.text(0.06,0.5, r"$\tilde{\nu}$ (cm$^{-1}$)", ha='center', va='center',rotation='vertical')

#fig3, axarr3 = plt.subplots(figRows,figCols,sharex='col',sharey='row',figsize=(4.3,6))
#fig3.subplots_adjust(wspace=0.10,hspace=0.35,left=0.12,right=0.98,top=0.93,bottom=0.1)
#fig3.text(0.5,0.02, "Time (ns)",ha='center', va='center') 
#fig3.text(0.06,0.5, r"Shift due to %s (cm$^{-1}$)"%field.replace('_',' ') , ha='center', va='center',rotation='vertical')

equilTime= 00 
forceDir = "force_calc"
for index,molec in enumerate(molecList) :
    ax = axarr[index%figRows,index/figRows]

    srfFile = "SNase_%s/%s/SNase_%s.solvent_rxn_field.projected.xvg"%(molec,forceDir,molec)
    pcfFile = "SNase_%s/%s/SNase_%s.protein_field.projected.xvg"%(molec,forceDir,molec)

    if not os.path.isfile(srfFile) or not os.path.isfile(pcfFile) : 
        print "File not found for %s"%molec
        sys.exit() 

    #srfData = np.genfromtxt(srfFile)  ## srf := Solvent Reaction Field 
    #pcfData = np.genfromtxt(pcfFile)  ## pcf := Protein Coloumb Field 

    #srfData *= (2.57 * -0.67)         ## kbT/eA -> MV/cm & MV/cm -> cm^-1
    #pcfData *= (2.57 * -0.67) 
    
    #srfData = srfData[int(float(equilTime)/simTime *len(srfData)):]
    #pcfData = pcfData[int(float(equilTime)/simTime *len(pcfData)):]


    ######
    ###   Reconstruct external field with "weighted" components
    ######
    #c = 1.00   ##Free parameter, weight of relative contribution from pcf
    #refData = srfData + c*pcfData   ## ref := Reconstructed External Field


    ######
    ###   Bin data and fit polynomial/guassian to histogram 
    ######
    if not os.path.isdir("data_fitting") : 
        os.mkdir('data_fitting') 

    #file="data_fitting/%s_solvent_rxn_field"%(molec) 
    #np.savetxt(file+".xvg",srfData)
    #command = "~/normal_distribution/tiltAngle -f %s.xvg -o %s.out -p %s.poly -g %s.gaus -t 25 --overwrite > /dev/null 2>&1 "%(file,file,file,file)
    #os.system(command)

    #file="data_fitting/%s_protein_field"%(molec) 
    #np.savetxt(file+".xvg",pcfData)
    #command = "~/normal_distribution/tiltAngle -f %s.xvg -o %s.out -p %s.poly -g %s.gaus -t 25 --overwrite > /dev/null 2>&1 "%(file,file,file,file)
    #os.system(command)

    #file="data_fitting/%s_reconstructed_field"%(molec) 
    #np.savetxt(file+".xvg",refData)
    #command = "~/normal_distribution/tiltAngle -f %s.xvg -o %s.out -p %s.poly -g %s.gaus -t 25 --overwrite > /dev/null 2>&1 "%(file,file,file,file)
    #os.system(command)

    ######
    ###   Read in polynomial/guassian fits
    ######
    srfData = np.genfromtxt("data_fitting/%s_solvent_rxn_field.poly"%(molec)) 
    pcfData = np.genfromtxt("data_fitting/%s_protein_field.poly"%(molec)) 
    #refData = np.genfromtxt("data_fitting/%s_reconstructed_field.poly"%(molec)) 
    refData = srfData

    ######
    ###   Normalize polynomial/guassian fits so that max of REF is 1
    ###### 
    #refMax = np.max(refData[:,1]) 
    #refData /= refMax

    srfMax  = np.max(srfData[:,1]) 
    pcfMax  = np.max(pcfData[:,1]) 

    if srfMax > pcfMax : 
        srfData[:,1] /= (srfMax * 2) 
        pcfData[:,1] /= (srfMax * 2) 
    else : 
        pcfData[:,1] /= (pcfMax * 2) 
        srfData[:,1] /= (pcfMax * 2) 

    ######
    ###   Plot polynomial/gaussian fits
    ######
    ax.plot(srfData[:,0],srfData[:,1],color=colorDict[molec],linestyle=':') 
    ax.plot(pcfData[:,0],pcfData[:,1],color=colorDict[molec],linestyle='--') 
    #ax.plot(refData[:,0],refData[:,1],color=colorDict[molec]) 

    ######
    ###   Read in experimental spectra 
    ######
    expFile = 'Exp_data/%s.dat'%molec
    if not os.path.isfile(expFile) : 
        print "Experimental spectra not found for %s"%molec
        continue 
    expData = np.genfromtxt(expFile) 

    ######
    ###   Plot experimental spectra    
    ######
    ax2 = ax.twiny() 
    ax2.plot(expData[:,0],expData[:,1],color=colorDict[molec],linestyle='-')


    ######
    ###   Set plot parameters
    ######
    ax2.set_xlim([2145,2180]) 
    ax.set_xlim([-50,100]) 

fig.savefig('figures/both_fields.png',format='png') 








    

