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
#"I92X",
"A109X",
#"V104X",
"N118X"   
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

figRows,figCols = 4,3
simTime =100 ##ns

#field = 'solvent_rxn_field'
srfDict, pcfDict = {}, {} 
srfStdDict, pcfStdDict = {}, {} 

equilTime= 60 
for field in ['total_field','external_field','solvent_rxn_field','protein_field' ] : 
    index=0
    fieldAccum = [] 
    stdAccum = [] 
    absMaxAccum = [] 

    fig, ax = plt.subplots(1,1,figsize=(4.3,3) )
    fig.subplots_adjust(left=0.20,bottom=0.13,right=0.95,top=0.95)
    fig.text(0.5,0.04, r"Shift due to %s (cm$^{-1}$)"%field.replace('_',' ') , ha='center', va='center') 
    fig.text(0.08,0.5, r"Occurances", ha='center', va='center',rotation='vertical') 
    
    fig2, ax2 = plt.subplots(1,1,figsize=(4.3,3) )
    fig2.subplots_adjust(left=0.20,bottom=0.13,right=0.95,top=0.95)
    fig2.text(0.5,0.02, r"Shift due to %s (cm$^{-1}$)"%field.replace('_',' '),ha='center', va='center') 
    fig2.text(0.06,0.5, r"$\tilde{\nu}$ (cm$^{-1}$)", ha='center', va='center',rotation='vertical')

    fig3, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row' )
    fig3.subplots_adjust(wspace=0.10,hspace=0.35,left=0.12,right=0.98,top=0.93,bottom=0.1)
    fig3.text(0.5,0.02, "Time (ns)",ha='center', va='center') 
    fig3.text(0.06,0.5, r"Shift due to %s (cm$^{-1}$)"%field.replace('_',' ') , ha='center', va='center',rotation='vertical')


    for index,molec in enumerate(molecList) :
        ax3 = axarr[index/figCols,index%figCols]

        timeFile = "SNase_%s/force_calc/SNase_%s.%s.projected.xvg"%(molec,molec,field)

        data = np.genfromtxt(timeFile)
        data = data[int(float(equilTime)/simTime *len(data)):]

        file="data_fitting/%s_%s"%(molec,field)
        np.savetxt(file+".xvg",data)
        command = "~/normal_distribution/tiltAngle -f %s.xvg -o %s.out -p %s.poly -g %s.gaus -t 25 --overwrite"%(file,file,file,file)
        os.system(command)
        data = np.genfromtxt('%s.gaus'%file)

        datafile = "data_fitting/%s_%s.out"%(molec,field)
        datafile2= "data_fitting/%s_%s.gaus"%(molec,field)

        try : 
            data = np.genfromtxt(datafile) 
            data2 = np.genfromtxt(datafile2) 
            timeData = np.genfromtxt(timeFile) 
        except IOError : 
            print "Skipping %s"%molec
            continue
    
        #meanField = data2[np.argmax(data2[:,1]),0]
        meanField = np.average(timeData[int(float(equilTime)/simTime *len(timeData)):])
        stdField = np.std(timeData[int(float(equilTime)/simTime *len(timeData)):]) 
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
            srfStdDict[molec] = stdField
        elif field == 'protein_field' : 
            pcfDict[molec] = meanField 
            pcfStdDict[molec] = stdField
        fieldAccum.append(meanField) 
        stdAccum.append(stdField) 
        absMaxAccum.append(absMaxDict[molec][0]) 

        ax2.scatter(meanField,absMaxDict[molec][0],color=nameToColorKeys[molec] ) 
    
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

plt.close(fig)
plt.close(fig2)
plt.close(fig3)

fig,ax = plt.subplots(1,1) 
fig2,ax2 = plt.subplots(1,1)
srfAccum, pcfAccum, freqAccum = [], [],[] 
fwhmAccum,errAccum = [],[]
fudgeFactor = 0.25
#for fudgeFactor in np.linspace(0.2,0.40,11) : 
for molec in srfDict : 
    #print molec, srfDict[molec], pcfDict[molec] 

    err = np.sqrt(srfStdDict[molec]**2 + fudgeFactor*(pcfStdDict[molec])**2)
    #err = srfStdDict[molec]
    ax.scatter(srfDict[molec]+pcfDict[molec]*fudgeFactor,absMaxDict[molec][0],color=nameToColorKeys[molec],s=5.0) 
    ax.errorbar(srfDict[molec]+pcfDict[molec]*fudgeFactor,absMaxDict[molec][0],yerr=absMaxDict[molec][2],xerr=err,color=nameToColorKeys[molec],capsize=2) 
    ax.scatter(srfDict[molec],absMaxDict[molec][0],color=nameToColorKeys[molec],alpha=0.25,s=5.0) 

    ax2.scatter(err,absMaxDict[molec][1],color=nameToColorKeys[molec]) 

    srfAccum.append(srfDict[molec]) 
    pcfAccum.append(pcfDict[molec]) 
    freqAccum.append(absMaxDict[molec][0]) 
    fwhmAccum.append(absMaxDict[molec][1]) 
    errAccum.append(err) 

xs = np.array(srfAccum) + np.array(pcfAccum)*fudgeFactor
slope, intercept, r_value, p_value, std_error = linregress(xs,freqAccum) 
xs = np.linspace(np.min(xs),np.max(xs),100) 
print "%0.2f\t%0.3f"%(fudgeFactor,r_value)
ax.plot(xs,slope*xs+intercept,'k--',label="r = %0.3f"%r_value) 

slope, intercept, r_value, p_value, std_error = linregress(errAccum,fwhmAccum) 
xs = np.linspace(np.min(errAccum),np.max(errAccum),100) 
print "FWHM: \t%0.3f"%(r_value)
ax2.plot(xs,slope*xs+intercept,'k--',label="r = %0.3f"%r_value) 

ax.set_xlabel(r"Shift due to Field (cm$^{-1}$)") 
ax.set_ylabel(r"Experimental frequency (cm$^{-1}$)") 
fig.legend()

ax2.set_xlabel(r"Standard deviation of Field (cm$^{-1}$)") 
ax2.set_ylabel(r"Experimental FWHM (cm$^{-1}$)") 
fig2.legend()

fig.savefig('figures/pcf+srf.png',format='png',dpi=500) 
fig2.savefig('figures/fwhm_vs_std.png') 




