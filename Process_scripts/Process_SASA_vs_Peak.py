import numpy as np 
import matplotlib.pyplot as plt 
from scipy.stats import linregress
import os
from matplotlib import rc_file

ftlsdata = 'Exp_data/ftls.dat'
rcFile='rc_files/paper.rc'
rc_file(rcFile) 

molecList = {
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
        "N118X"
}

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

peakDict = {}
with open('Exp_data/abs_data.dat') as f :
    for line in f :
        if not line.startswith('#') :
            key = line.split()[0]
            peak, FWHM, error = line.split()[1:4]
            peak, FWHM, error = float(peak), float(FWHM), float(error)
            peakDict[key] = [peak,FWHM,error]

for sasa in ['cnc_area','sidechain_area','thiocyanate','nitrile_area'] : 
    peakAccum, sasaAccum = [],[]
    for molec in molecList : 
        peak, error = peakDict[molec][0], peakDict[molec][2]
    
        sasaFile = "SNase_%s/sasa/%s.xvg"%(molec,sasa) 
    
        if not os.path.isfile(sasaFile) : 
            print "%s %s not found."%(sasaFile,molec) 
            continue 
    
        headlines = 0 
        with open(sasaFile) as f : 
            for line in f.readlines() : 
                if line.startswith('#') or line.startswith('@') : 
                    headlines += 1 
                else : 
                    break 
        try : 
            data = np.genfromtxt(sasaFile,skip_header=headlines) 
        except ValueError : 
            print "Error import %s %s" %(molec, sasa) 
            continue 
        
        sasaAvg = np.average(data[:,2]) 
    
        plt.scatter(sasaAvg,peak,color=nameToColorKeys[molec]) 
    
        peakAccum.append(peak)
        sasaAccum.append(sasaAvg) 
    
    
    slope, intercept, r_value, p_value, std_error = linregress(sasaAccum,peakAccum) 
    print "r = %0.3f"%r_value
    
    x = np.linspace(min(sasaAccum),max(sasaAccum),100) 
    plt.plot(x,slope*x+intercept,'k--',label='%0.3f'%r_value) 
    
    plt.xlabel("SASA (%s) (nm$^2$)"%(sasa.split('_')[0]) ) 
    plt.ylabel("Experimental Frequency (cm$^{-1}$)") 
    plt.legend()
    plt.savefig('figures/peak_vs_sasa_%s.png'%(sasa.split('_')[0]),format='png',dpi=500) 
    plt.close()

        


