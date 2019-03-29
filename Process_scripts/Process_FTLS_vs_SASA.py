import numpy as np 
import matplotlib.pyplot as plt 
from scipy.stats import linregress
import os
from matplotlib import rc_file

ftlsdata = 'Exp_data/ftls.dat'
rcFile='rc_files/paper.rc'
rc_file(rcFile) 

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

ftlsDict = {}
with open(ftlsdata) as f : 
    for line in f.readlines() : 
        key,value, error = line.split()
        if value == "nan" : continue 
        ftlsDict[key] = [float(value), float(error)] 

for sasa in ['cnc_area','sidechain_area','thiocyanate','nitrile_area'] : 
    ftlsAccum, sasaAccum = [],[]
    for molec in ftlsDict : 
        ftls, error = ftlsDict[molec]
    
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
    
        plt.scatter(sasaAvg,ftls,color=nameToColorKeys[molec]) 
    
        ftlsAccum.append(ftls)
        sasaAccum.append(sasaAvg) 
    
    
    slope, intercept, r_value, p_value, std_error = linregress(sasaAccum,ftlsAccum) 
    print "r = %0.3f"%r_value
    
    x = np.linspace(min(sasaAccum),max(sasaAccum),100) 
    plt.plot(x,slope*x+intercept,'k--',label='%0.3f'%r_value) 
    
    plt.xlabel("SASA (%s) (nm$^2$)"%(sasa.split('_')[0]) ) 
    plt.ylabel("FTLS") 
    plt.legend()
    plt.savefig('figures/ftls_vs_sasa_%s.png'%(sasa.split('_')[0]),format='png',dpi=500) 
    plt.close()

        


