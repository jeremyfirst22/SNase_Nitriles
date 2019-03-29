import glob 
import numpy as np 
import matplotlib.pyplot as plt 
import os
from os import sys
from matplotlib.colors import LogNorm
import matplotlib.lines as mlines 
from matplotlib import rc_file
from scipy.stats import linregress

figCols=1
figRows=1

nameToColorKeys = {
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
molecList = [
"V23X",
"L25X",
"L38X",
"A58X",
"T62X",
"V66X",
"A90X",
"I92X",
"A109X",
#"V104X",
#"N118X"   
]

absMax = {}  
with open('Exp_data/sasa_abs_data.dat') as f : 
    for line in f : 
        if not line.startswith('#') : 
            key = line.split()[0] 
            value = float(line.split()[2]) 
            absMax[key] = value 


rcFile='rc_files/paper.rc'
rc_file(rcFile) 

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

figRows,figCols = 4,3
simTime =100 ##ns

#field = 'solvent_rxn_field'
srfDict, pcfDict = {}, {} 
for field in ['total_field','external_field','solvent_rxn_field','protein_field' ] : 
    index=0
    fieldAccum = [] 
    absMaxAccum = [] 

    fig, ax = plt.subplots(1,1,figsize=(4.3,3) )
    fig.subplots_adjust(left=0.20,bottom=0.13,right=0.95,top=0.95)
    fig.text(0.5,0.04, "Field ()", ha='center', va='center') 
    fig.text(0.08,0.5, r"Occurances", ha='center', va='center',rotation='vertical') 
    
    fig2, ax2 = plt.subplots(1,1,figsize=(4.3,3) )
    fig2.subplots_adjust(left=0.20,bottom=0.13,right=0.95,top=0.95)
    fig2.text(0.5,0.02, "Shift due to Field (cm$^{-1}$)",ha='center', va='center') 
    fig2.text(0.06,0.5, r"$\tilde{\nu}$ (cm$^{-1}$)", ha='center', va='center',rotation='vertical')

    fig3, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row' )
    fig3.subplots_adjust(wspace=0.10,hspace=0.35,left=0.12,right=0.98,top=0.93,bottom=0.1)
    fig3.text(0.5,0.02, "Time (ns)",ha='center', va='center') 
    fig3.text(0.06,0.5, r"Shift due to field (cm$^{-1}$)", ha='center', va='center',rotation='vertical')


    for index,molec in enumerate(molecList) :
        ax3 = axarr[index/figCols,index%figCols]

        datafile = "SNase_%s/force_calc/%s.out"%(molec,field) 
        datafile2= "SNase_%s/force_calc/%s.poly"%(molec,field)

        timeFile = "SNase_%s/force_calc/SNase_%s.%s.projected.xvg"%(molec,molec,field)
    
        try : 
            data = np.genfromtxt(datafile) 
            data2 = np.genfromtxt(datafile2) 
            timeData = np.genfromtxt(timeFile) 
        except IOError : 
            print "Skipping %s"%molec
            continue
    
        #meanField = data2[np.argmax(data2[:,1]),0]
        meanField = np.average(timeData) 
        meanField *= 2.57 ##kbT/eA -> MV/cm
        meanField *= -0.67 ## MV/cm -> cm^-1

        data[:,0] *= (2.57 * -0.67) 
        data2[:,0] *= (2.57 * -0.67) 
        timeData *= (2.57 * -0.67)

        time = np.linspace(0,simTime,len(timeData) ) 
        ax3.scatter(time,timeData,s=0.05,color=nameToColorKeys[molec],marker='.')  
        ax3.set_title(molec,color=nameToColorKeys[molec]) 

        max = np.max(data2[:,1])
    #    data[:,1] /= max
    #    data2[:,1] /=max
    
        ax.scatter(data[:,0],data[:,1],s=0.1,color=nameToColorKeys[molec]) 
        ax.plot(data2[:,0],data2[:,1],color=nameToColorKeys[molec]) 
        ax.axvline(np.average(timeData),color=nameToColorKeys[molec],linestyle='--') 
    
        print molec, meanField
        if field == 'solvent_rxn_field' : 
            srfDict[molec] = meanField 
        elif field == 'protein_field' : 
            pcfDict[molec] = meanField 

        ax2.scatter(meanField,absMax[molec],color=nameToColorKeys[molec] ) 
    
        fieldAccum.append(meanField) 
        absMaxAccum.append(absMax[molec]) 
        index +=1
    
    slope,intercept,r_value,p_value,std_error = linregress(fieldAccum, absMaxAccum)
    print "r = %f, p = %f"%(r_value,p_value)
    
    xs = np.linspace(np.min(fieldAccum), np.max(fieldAccum),100 )
    ys = slope * xs + intercept
    ax2.plot(xs, ys, label="r = %0.3f"%r_value,color='k')
    ax2.legend(loc=2) 
    
    fig.savefig('figures/%s.png'%field,format='png') 
    fig2.savefig('figures/%s_vs_peak.png'%field,format='png') 
    fig3.savefig('figures/time_vs_%s.png'%field,format='png') #,dpi=1500) 

fig,ax = plt.subplots(1,1) 
srfAccum, pcfAccum, freqAccum = [], [],[] 
for molec in srfDict : 
    print molec, srfDict[molec], pcfDict[molec] 

    ax.scatter(srfDict[molec]+pcfDict[molec]*0.2,absMax[molec],color=nameToColorKeys[molec]) 
    ax.scatter(srfDict[molec],absMax[molec],color=nameToColorKeys[molec],alpha=0.25) 

    srfAccum.append(srfDict[molec]) 
    pcfAccum.append(pcfDict[molec]) 
    freqAccum.append(absMax[molec]) 

xs = np.array(srfAccum) + np.array(pcfAccum)*0.2
slope, intercept, r_value, p_value, std_error = linregress(xs,freqAccum) 
xs = np.linspace(np.min(xs),np.max(xs),100) 
ax.plot(xs,slope*xs+intercept,'k--',label="r = %0.3f"%r_value) 
fig.legend()

fig.savefig('figures/pcf+srf.png',format='png',dpi=500) 




