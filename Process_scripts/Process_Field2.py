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

absdata = 'Exp_data/abs_data2.dat'

peakDict, fwhmDict = {}, {}
with open(absdata) as f :
    for line in f.readlines() :
        if line.startswith('#') : continue
        key,peak,peakError,fwhm,fwhmError = line.split()
        if fwhm == "nan" : continue
#        value = float(value) * -1
        fwhmDict[key] = [float(fwhm),float(fwhmError) ]
        peakDict[key] = [float(peak),float(peakError) ]

rcFile='rc_files/paper.rc'
rc_file(rcFile) 

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

figRows,figCols = 5,2
simTime =100 ##ns

#field = 'solvent_rxn_field'
srfDict, pcfDict = {}, {} 
srfStdDict, pcfStdDict = {}, {} 


equilTime= 0 
#for field in ['solvent_rxn_field','protein_field','external_field'] : 
for field in ['solvent_rxn_field','protein_field'] : #,'external_field'] : 
    index=0
    fieldAccum = [] 
    stdAccum = [] 
    peakAccum = [] 
    fwhmAccum = [] 

    left, right = 0.15, 0.95
    bottom, top = 0.13, 0.98
    wspace, hspace = 0.1, 0.0

    fig, axarr = plt.subplots(figRows,figCols,sharex='all',sharey='all',figsize=(3,3) )
    fig.subplots_adjust(left=left,bottom=bottom,right=right,top=top,hspace=hspace,wspace=wspace)
    fig.text((right-left)/2+left,0.03, r"%s ($\frac{k_BT}{e^- \rm{\AA}}$)"%field.capitalize().replace('_',' ')  , ha='center', va='center') 
    fig.text(0.04,(top-bottom)/2+bottom, r"Probability (a.u.)", ha='center', va='center',rotation='vertical') 
    axarr[0,0].invert_xaxis()

    left, right = 0.18, 0.95
    bottom, top = 0.18, 0.95
    
    fig2, ax2 = plt.subplots(1,1,figsize=(3,2) )
    fig2.subplots_adjust(left=left,bottom=bottom,right=right,top=top)
    fig2.text((right-left)/2+left,0.01, r"%s (cm$^{-1}$)"%field.capitalize().replace('_',' '),ha='center', va='bottom') 
    fig2.text(0.01,(top-bottom)/2+bottom, r"$\tilde{\nu}$ (cm$^{-1}$)", ha='left', va='center',rotation='vertical')
    ax2.invert_xaxis() 

    fig4, ax4 = plt.subplots(1,1,figsize=(3,2) )
    fig4.subplots_adjust(left=left,bottom=bottom,right=right,top=top)
    fig4.text((right-left)/2+left,0.01, r"Std. dev. of %s (cm$^{-1}$)"%field.capitalize().replace('_',' '),ha='center', va='bottom') 
    fig4.text(0.01,(top-bottom)/2+bottom, r"FWHM (cm$^{-1}$)", ha='left', va='center',rotation='vertical')

    left, right = 0.13, 0.95
    bottom, top = 0.13, 0.98
    wspace, hspace = 0.1, 0.0

    fig3, axarr3 = plt.subplots(figRows,figCols,sharex='all',sharey='all',figsize=(3,3))
    fig3.subplots_adjust(wspace=wspace,hspace=hspace,left=left,right=right,top=top,bottom=bottom)
    fig3.text((right-left)/2+left,0.02, "Time (ns)",ha='center', va='center') 
    fig3.text(.06,(top-bottom)/2+bottom, r"%s ($\frac{k_BT}{e^- \rm{\AA}}$)"%field.replace('_',' ') , ha='center', va='center',rotation='vertical')
    axarr3[0,0].invert_yaxis() 


    for index,molec in enumerate(molecList) :
        ax = axarr[index%figRows,index/figRows]
        ax3 = axarr3[index%figRows,index/figRows]

        timeFile = "SNase_%s/force_calc_ca/SNase_%s.%s.projected.xvg"%(molec,molec,field)

        if not os.path.isfile(timeFile) : 
            print "%s not found. Skipping"%timeFile
            continue 
        
        data = np.genfromtxt(timeFile)
        data = data[int(float(equilTime)/simTime *len(data)):]

        file="data_fitting/%s_%s"%(molec,field)
        np.savetxt(file+".xvg",data)
        command = "~/normal_distribution/tiltAngle -f %s.xvg -o %s.out -p %s.poly -g %s.gaus -t 25 --overwrite >/dev/null 2>&1"%(file,file,file,file)
        os.system(command)

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
        meanField = np.average(timeData[int(float(equilTime)/simTime *len(timeData)):])
        stdField = np.std(timeData[int(float(equilTime)/simTime *len(timeData)):]) 
        #meanField *= 2.57 ##kbT/eA -> MV/cm
        #meanField *= -0.67 ## MV/cm -> cm^-1

        #data[:,0] *= (2.57 * -0.67) 
        #data2[:,0] *= (2.57 * -0.67) 
        #timeData *= (2.57 * -0.67)

        time = np.linspace(0,simTime,len(timeData) ) 
#        ax3.scatter(time,timeData,s=0.05,color=colorDict[molec],marker='.')  
#        ax3.text(0.95,0.95,"%s"%molec,va='top',ha='right',color=colorDict[molec],transform=ax3.transAxes)

        #ax3.set_ylim([-60,60]) 

        max = np.max(data2[:,1])
        data[:,1] /= max
        data2[:,1] /=max
    
        #ax.scatter(data[:,0],data[:,1],s=0.1,color=colorDict[molec]) 
        binwidth = data[0,0] - data[1,0] 
#        ax.plot(data2[:,0],data2[:,1],color=colorDict[molec],linewidth=0.5 )
#        ax.bar(data[::1,0],data[::1,1],width=binwidth,color='gray',alpha=0.5 )
#        ax.text(0.95,0.85,"%s"%molec,va='top',ha='right',color=colorDict[molec],transform=ax.transAxes)
        #ax.axvline(np.average(timeData),color=colorDict[molec],linestyle='--',linewidth=1) 
        #ax.set_xlim([50,-50]) 
    
        print molec, meanField

        if field == 'solvent_rxn_field' : 
            srfDict[molec] = meanField 
            srfStdDict[molec] = stdField
        elif field == 'protein_field' : 
            pcfDict[molec] = meanField 
            pcfStdDict[molec] = stdField
        fieldAccum.append(meanField) 
        stdAccum.append(stdField) 
        peakAccum.append(peakDict[molec][0]) 
        fwhmAccum.append(fwhmDict[molec][0]) 

#        ax2.scatter(meanField,absMaxDict[molec][0],color=colorDict[molec]) 
#        ax4.scatter(stdField,absMaxDict[molec][1],color=colorDict[molec]) 
    
        index +=1
    
#    slope,intercept,r_value,p_value,std_error = linregress(fieldAccum, peakAccum)
#    print "r = %f, p = %f"%(r_value,p_value)
#    xs = np.linspace(np.min(fieldAccum), np.max(fieldAccum),100 )
#    ys = slope * xs + intercept
#    ax2.plot(xs, ys,linestyle='--',color='k')
#    ax2.text(0.1,0.9,"r = %0.2f"%r_value,transform=ax2.transAxes,va='top',ha='left')
#
#    slope,intercept,r_value,p_value,std_error = linregress(stdAccum, fwhmAccum)
#    print "r = %f, p = %f"%(r_value,p_value)
#    xs = np.linspace(np.min(stdAccum), np.max(stdAccum),100 )
#    ys = slope * xs + intercept
#    ax4.plot(xs, ys,linestyle='--',color='k')
#    ax4.text(0.1,0.9,"r = %0.2f"%r_value,transform=ax4.transAxes,va='top',ha='left')
#    #ax2.legend(loc=2) 
#
#    ax.set_yticks([])
#
#    
#    fig2.savefig('figures/%s_vs_peak.png'%field,format='png',dpi=500) 
#    fig4.savefig('figures/std_%s_vs_fwhm.png'%field,format='png',dpi=500) 
#    fig3.savefig('figures/time_vs_%s.png'%field,format='png',dpi=1500) 
#    fig.savefig('figures/%s.png'%field,format='png',dpi=500) 

plt.close(fig)
plt.close(fig2)
plt.close(fig3)

left, right = 0.15, 0.80
bottom, top = 0.18, 0.98

fig, ax = plt.subplots(1,1,figsize=(3.42,2.25)) 
fig.subplots_adjust(left=left,bottom=bottom,right=right,top=top,hspace=hspace,wspace=wspace)
fig.text((right-left)/2+left,0.03, r"Avg. of PCF + SRF ($\frac{k_BT}{e^- \rm{\AA}}$)", ha='center', va='center') 
fig.text(0.04,(top-bottom)/2+bottom, r"$\nu$ (cm$^{-1}$)", ha='center', va='center',rotation='vertical') 
ax.invert_xaxis() 

fig2, ax2= plt.subplots(1,1,figsize=(3.42,2.25)) 
fig2.subplots_adjust(left=left,bottom=bottom,right=right,top=top,hspace=hspace,wspace=wspace)
fig2.text((right-left)/2+left,0.01, r"$(\sigma_{\text{SRF}}^2 + \sigma_{\text{PCF}}^2)^{\frac{1}{2}} \left ( \frac{k_BT}{e^- \rm{\AA}} \right ) $", ha='center', va='bottom') 
fig2.text(0.01,(top-bottom)/2+bottom, r"FWHM (cm$^{-1}$)", ha='left', va='center',rotation='vertical') 

srfAccum, pcfAccum, freqAccum = [], [],[] 
fwhmAccum,stdAccum = [],[]
#for fudgeFactor in np.linspace(0.0,1.50,11) : 
for fudgeFactor in [1.00] : 
    for molec in srfDict : 
        #print molec, srfDict[molec], pcfDict[molec] 
    
        std = np.sqrt(srfStdDict[molec]**2 + fudgeFactor*(pcfStdDict[molec])**2)
        ax.scatter(srfDict[molec]+pcfDict[molec]*fudgeFactor,peakDict[molec][0],color=colorDict[molec]) 
#        ax.errorbar(srfDict[molec]+pcfDict[molec]*fudgeFactor,absMaxDict[molec][0],yerr=absMaxDict[molec][2],xerr=err,color=colorDict[molec],capsize=2) 
        #ax.scatter(srfDict[molec],absMaxDict[molec][0],color=colorDict[molec],alpha=0.25,s=5.0) 
    
        ax2.scatter(std,fwhmDict[molec][0],color=colorDict[molec],label=molec) 
        ax2.errorbar(std,fwhmDict[molec][0],yerr=fwhmDict[molec][1],color=colorDict[molec])
    
        srfAccum.append(srfDict[molec]) 
        pcfAccum.append(pcfDict[molec]) 
        freqAccum.append(peakDict[molec][0]) 
        fwhmAccum.append(fwhmDict[molec][0]) 
        stdAccum.append(std) 
    
xs = np.array(srfAccum) + np.array(pcfAccum)*fudgeFactor
slope, intercept, r_value, p_value, std_error = linregress(xs,freqAccum) 
xs = np.linspace(np.min(xs),np.max(xs),100) 
print "%0.2f\t%0.3f"%(fudgeFactor,r_value)
ax.plot(xs,slope*xs+intercept,'k--',) 
ax.text(0.1,0.9,"$r$ = %0.2f"%r_value,transform=ax.transAxes,ha='left',va='top')

slope, intercept, r_value, p_value, std_error = linregress(stdAccum,fwhmAccum) 
xs = np.linspace(np.min(stdAccum),np.max(stdAccum),100) 
print "FWHM: \t%0.3f"%(r_value)
ax2.plot(xs,slope*xs+intercept,'k--') 
ax2.text(0.05,0.95,"$r$ = %0.2f"%r_value,transform=ax2.transAxes,ha='left',va='top')
#fig2.text(0.96,0.27,r"$F_{SCF}+$%1.2f$\cdot F_{PCF}$"%(fudgeFactor),ha='right',va='center') 
#ax2.legend(loc=4)
ax2.legend(loc=(0.97,0.12))

fig.savefig('figures/pcf+srf_vs_peak.png',format='png',dpi=500) 
fig2.savefig('figures/fwhm_vs_std.png',dpi=500) 




