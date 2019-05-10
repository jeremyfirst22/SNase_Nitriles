import glob 
import numpy as np 
import matplotlib.pyplot as plt 
import os
from os import sys
from matplotlib.colors import LogNorm
import matplotlib.lines as mlines 
from matplotlib import rc_file
from scipy.stats import linregress
from matplotlib import collections as matcoll

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

ftlsDict = {}

ftlsdata = 'Exp_data/ftls_fits.dat'
with open(ftlsdata) as f :
    for line in f.readlines() :
        key,value, error = line.split()
        if value == "nan" : continue
        #value = float(value) * -1
        ftlsDict[key] = [float(value), float(error)]


rcFile='rc_files/paper.rc'
rc_file(rcFile) 

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

figRows,figCols = 5,2
simTime =100 ##ns

#field = 'solvent_rxn_field'
srfDict, pcfDict = {}, {} 
srfStdDict, pcfStdDict = {}, {} 

def gauss (a,b,c,x) : 
    return a*np.exp(-(x-b)**2/(c**2))
def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return (average, np.sqrt(variance))

equilTime= 0
#xmin,xmax, xstep = -60,100,500

def find_nearest(array, value) : 
    array = np.asarray(array)
    idx = np.argmin(np.abs(array - value))
    return idx

def full_width_half_max(x, y) : 
    peakIndex = np.argmax(y)
    minIndex = find_nearest(y[:peakIndex],y[peakIndex] * 0.5) 
    maxIndex = find_nearest(y[peakIndex:],y[peakIndex] * 0.5) + peakIndex

    return np.abs(x[maxIndex] - x[minIndex]) #width at half peak


for field in ['external_field'] : 
    index=0
    fieldAccum = [] 
    stdAccum = [] 
    fwhmAccum = [] 
    absMaxAccum = [] 

    left, right = 0.08, 0.99
    bottom, top = 0.11, 0.99
    wspace, hspace = 0.05, 0.0

    fig, axarr = plt.subplots(figRows,figCols,sharex='all',sharey='all',figsize=(3.42,3) )
    fig.subplots_adjust(left=left,bottom=bottom,right=right,top=top,hspace=hspace,wspace=wspace)
    fig.text((right-left)/2+left,0.01, r"Wavenumbers (cm$^{-1}$)", ha='center', va='bottom') 
    fig.text(0.01,(top-bottom)/2+bottom, r"Absorbance or probability (a.u.)", ha='left', va='center',rotation='vertical') 

    left, right = 0.13, 0.95 
    bottom, top = 0.13, 0.98 
    hspace, wspace = 0, 0.1

    fig2, ax2 = plt.subplots(1,1,figsize=(3.42,2) )
    fig2.subplots_adjust(left=left,bottom=bottom,right=right,top=top,hspace=hspace,wspace=wspace)
    fig2.text((right-left)/2+left,0.03, r"Std. Dev. of broadened field (cm$^{-1}$)", ha='center', va='center') 
    fig2.text(0.04,(top-bottom)/2+bottom, r"FWHM (cm$^{-1}$)", ha='center', va='center',rotation='vertical') 

    for index,molec in enumerate(molecList) :
        ax = axarr[index%figRows,index/figRows]
#        ax3 = axarr3[index%figRows,index/figRows]

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
        except IOError : 
            print "Skipping %s"%molec
            continue

        expFile = 'Exp_data/%s.dat'%molec
        if not os.path.isfile(expFile) :
            print "Experimental spectra not found for %s"%molec
            continue
        expData = np.genfromtxt(expFile)
        ax.plot(expData[:,0],expData[:,1],color='k',linestyle='--',linewidth=1) 

        data2[:,1] /= np.max(data2[:,1]) 

        avg, std = weighted_avg_and_std(data[:,0],data[:,1])
        #### Move average of Forces to zero
        data[:,0] -= data2[np.argmax(data2[:,1]),0] 
        data2[:,0] -= data2[np.argmax(data2[:,1]),0] 
        ###  Match standard deviations using FWHM of V66X as reference (3.375 cm^-1)
        #conversion = 1/12.50  ## original = 2.33 
        conversion = -0.3953 ## wavenumbers per kb T / e A
        data[:,0] *= conversion 
        data2[:,0] *= conversion 
        #data[:,0]  *= (-3.375)*conversion  
        #data2[:,0] *= (-3.375)*conversion  
        ### Normalize histogram height to 1
        data[:,1] /= np.max(data[:,1]) 

        ### Use FTLS to scale bandwidth  
        #conversion =  155.0 
        #minimumWidening = 0.90  
        #conversion = 0 
        #y1, y2 = 1.7229721945322476, 6.183000963739933 ## Standard deviation of V66 and N118
        #x1, x2 = ftlsDict["V66X"][0], ftlsDict["N118X"][0]   ## FTLS of V66 and N118
        #conversion = ( y2 - y1 ) / (x2 - x1 ) 
        #minimumWidening = y1 - conversion * x1
        
        conversion = 144 #190 
        minimumWidening = 1.71 - (conversion * -ftlsDict["V66X"][0])  ##1.71 matches V66 FWHM best. 

        bandwidth = conversion * -ftlsDict[molec][0] + minimumWidening

        print "\t%10s\t%10f\t%10f\t%10f\t%10f\t%10f"%(molec, bandwidth, ftlsDict[molec][0], conversion, minimumWidening, ftlsDict[molec][0] * conversion + minimumWidening)
        #bandwidth = -conversion * ftlsDict[molec][0]
        
        ### Broaden histogram bins with gaussian
        xmin, xmax, xstep = - 40,  40,500
        xs = np.linspace(xmin, xmax, xstep, dtype='float') 
        bindata = np.zeros_like(xs) 
        for trans in data : 
            a,b,c = trans[1],trans[0],bandwidth
            spec = gauss(a,b,c,xs)
            
            for k,x in enumerate(xs) : 
                bindata[k] += spec[k]

        bindata /= np.max(bindata)      ##Normalize to 1 

        max = np.max(data2[:,1])
        data[:,1] /= max
        data2[:,1] /=max

        #data[:,0] += data2[np.argmax(expData[:,1]), 0]
        #data2[:,0] += expData[np.argmax(expData[:,1]),0]#absMaxDict[molec][0]

        ### Move Peak of broadened data to peak of experiment
        xs += expData[np.argmax(expData[:,1]),0] - xs[np.argmax(bindata)]#absMaxDict[molec][0]
        data[:,0] += expData[np.argmax(expData[:,1]),0]#absMaxDict[molec][0]
        data2[:,0] += expData[np.argmax(expData[:,1]),0]#absMaxDict[molec][0]

        assert xs[np.argmax(bindata) ] == expData[np.argmax(expData[:,1]),0]


        ### Plot histogram
        binwidth = data[0,0] - data[1,0]
        ax.bar(data[::1,0],data[::1,1],width=binwidth,color='gray',alpha=0.5 )
        ax.plot(data2[:,0], data2[:,1],'k--',linewidth=0.25,zorder=10) 
    
        ax.plot(xs,bindata,color=colorDict[molec],zorder=1 )

    
        index +=1
        ax.set_yticks([]) 
        #ax.set_xticks([2145,2160,2175]) 
        ax.set_xlim([2140,2180]) 
        ax.set_xticks([2145,2160,2175])


        avgB, stdB= weighted_avg_and_std(xs, bindata) 

        _, stdExp = weighted_avg_and_std(expData[:,0], expData[:,1]) 

        fwhm = fwhmDict[molec][0]
        ax2.scatter(stdB, fwhm,color=colorDict[molec]) 
        #ax2.scatter(ftlsDict[molec][0], fwhm,color=colorDict[molec]) 

        #stdAccum.append(ftlsDict[molec][0]) 
        stdAccum.append(stdB) 
        label = molec 
        if molec in ["V66X","N118X"] : 
            print molec, full_width_half_max(xs, bindata) , fwhmDict[molec][0]
            label += "*"
        ax.text(0.01,0.88,"%s"%label,va='top',ha='left',color=colorDict[molec],transform=ax.transAxes)
        
        fwhmAccum.append(fwhm) 

    slope,intercept,r_value,p_value,std_error = linregress(stdAccum, fwhmAccum)
    print "r = %f, p = %f"%(r_value,p_value)
    xs = np.linspace(np.min(stdAccum), np.max(stdAccum),100 )
    ys = slope * xs + intercept
    ax2.plot(xs, ys,linestyle='--',color='k')
    ax2.text(0.1,0.9,"r = %0.2f"%r_value,transform=ax2.transAxes,va='top',ha='left')
    
    fig.savefig('figures/convoluted_%s.png'%field,format='png',dpi=500) 
    fig2.savefig('figures/convoluted_std_vs_fwhm.png',format='png',dpi=500) 



plt.close(fig)
#plt.close(fig2)
#plt.close(fig3)

sys.exit() 

