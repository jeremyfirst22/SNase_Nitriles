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

figCols=2
figRows=5

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

srfDict, pcfDict = {}, {} 
srfStdDict, pcfStdDict = {}, {} 

fig = plt.figure(figsize=(6.5,8))

left, right    = 0.08, 0.95 
bottom, top    = 0.06, 0.97 
hspace, wspace = 0.15, 0.20

plotArray = gridspec.GridSpec(2,2,wspace=wspace,hspace=hspace,left=left,right=right,bottom=bottom,top=top)

fig.text(0.01,top-(top-bottom-hspace/2)/4   , r"SRF $(\frac{k_B T}{e^{-}\AA})$",ha='left',va='center',rotation='vertical')
fig.text(0.01,bottom+(top-bottom-hspace/2)/4, r"PCF $(\frac{k_B T}{e^{-}\AA})$",ha='left',va='center',rotation='vertical')
fig.text(0.53,top-(top-bottom-hspace/2)/4   , r"Probability (a.u.)             ",ha='left',va='center',rotation='vertical')
fig.text(0.53,bottom+(top-bottom-hspace/2)/4, r"Probability (a.u.)            ",ha='left',va='center',rotation='vertical')

fig.text(right-(right-left-wspace/2)/4,0.50 , r"SRF $(\frac{k_B T}{e^{-}\AA})$",ha='center',va='bottom')
fig.text(right-(right-left-wspace/2)/4,0.01 , r"PCF $(\frac{k_B T}{e^{-}\AA})$",ha='center',va='bottom')
fig.text(left +(right-left-wspace/2)/4,0.015, r"Time (ns)                     ",ha='center',va='bottom')
fig.text(left +(right-left-wspace/2)/4,0.505, r"Time (ns)                     ",ha='center',va='bottom') 

fig.text(0.01,0.98,r"\textsf{A}",va='top',ha='left',fontsize=12)
fig.text(0.52,0.98,r"\textsf{B}",va='top',ha='left',fontsize=12)
fig.text(0.01,0.49,r"\textsf{C}",va='top',ha='left',fontsize=12)
fig.text(0.52,0.49,r"\textsf{D}",va='top',ha='left',fontsize=12)

for plotRow,field in enumerate(['solvent_rxn_field','protein_field']) : 
    pcfstdAccum = [] 
    srfstdAccum = [] 

    innerTime = gridspec.GridSpecFromSubplotSpec(figRows,figCols,subplot_spec=plotArray[plotRow,0],wspace=0.1,hspace=0.00)
    innerHist = gridspec.GridSpecFromSubplotSpec(figRows,figCols,subplot_spec=plotArray[plotRow,1],wspace=0.1,hspace=0.00)
    for index,molec in enumerate(molecList) :
        axTime = plt.Subplot(fig,innerTime[index%figRows,index/figRows])
        fig.add_subplot(axTime,sharex='all',sharey='all')
        axTime.text(0.96,0.93,"%s"%molec,va='top',ha='right'\
                ,color=colorDict[molec],transform=axTime.transAxes)

        axHist = plt.Subplot(fig,innerHist[index%figRows,index/figRows])
        fig.add_subplot(axHist,sharex='all',sharey='all')
        axHist.text(0.96,0.93,"%s"%molec,va='top',ha='right'\
                ,color=colorDict[molec],transform=axHist.transAxes)

        if not index%figRows + 1 == figRows :
            axTime.tick_params(axis="x",labelbottom=False)
            axHist.tick_params(axis="x",labelbottom=False)
        if index/figRows + 1 == figCols :
            axTime.tick_params(axis="y",labelleft=False)
        axHist.tick_params(axis="y",labelleft=False,left=False)

        axTime.set_ylim([25,-50]) 
        axHist.set_xlim([25,-50]) 

        axTime.set_yticks(np.arange(-40,30,20))
        axHist.set_xticks(np.arange(-40,30,20))

#        continue 

        timeFile = "SNase_%s/force_calc_ca/SNase_%s.%s.projected.xvg"%(molec,molec,field)

        if not os.path.isfile(timeFile) : 
            print "%s not found. Skipping"%timeFile
            continue 
        
        #data = np.genfromtxt(timeFile)
        #data = data[int(float(equilTime)/simTime *len(data)):]

#        file="data_fitting/%s_%s"%(molec,field)
#        np.savetxt(file+".xvg",data)
#        command = "~/normal_distribution/tiltAngle -f %s.xvg -o %s.out -p %s.poly -g %s.gaus -t 25 --overwrite >/dev/null 2>&1"%(file,file,file,file)
#        os.system(command)

        datafile = "data_fitting/%s_%s.out"%(molec,field)
        datafile2= "data_fitting/%s_%s.poly"%(molec,field)

        try : 
            data = np.genfromtxt(datafile) 
            data2 = np.genfromtxt(datafile2) 
            timeData = np.genfromtxt(timeFile) 
        except IOError : 
            print "Skipping %s"%molec
            continue
    
        #meanField = data2[np.argmax(data2[:,1]),0]
        #meanField *= 2.57 ##kbT/eA -> MV/cm
        #meanField *= -0.67 ## MV/cm -> cm^-1

        #data[:,0] *= (2.57 * -0.67) 
        #data2[:,0] *= (2.57 * -0.67) 
        #timeData *= (2.57 * -0.67)

        time = np.linspace(0,simTime,len(timeData) ) 
        axTime.scatter(time,timeData,s=0.05,color=colorDict[molec],marker='.')  

        #ax3.set_ylim([-60,60]) 

        max = np.max(data2[:,1])
        data[:,1] /= max
        data2[:,1] /=max
    
        #ax.scatter(data[:,0],data[:,1],s=0.1,color=colorDict[molec]) 
        binwidth = data[0,0] - data[1,0] 
        axHist.plot(data2[:,0],data2[:,1],color=colorDict[molec],linewidth=0.5 )
        axHist.bar(data[::1,0],data[::1,1],width=binwidth,color='gray',alpha=0.5 )
        #axHist.text(0.95,0.85,"%s"%molec,va='top',ha='right',color=colorDict[molec],transform=axHist.transAxes)
        axHist.axvline(np.average(timeData),color=colorDict[molec],linestyle='--',linewidth=1) 
        #ax.set_xlim([50,-50]) 

fig.savefig('figures/combined_figure_fields.png',format='png') 
plt.close(fig)

