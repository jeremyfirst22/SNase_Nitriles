import numpy as np 
import matplotlib.pyplot as plt 
import math
from scipy.stats import linregress 

file = 'data.txt'
outfile = '../ftls.dat'

with open(file) as f : 
    legend = f.readlines()[0]
    legend = legend.split() 

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
        "N118X":'#7E2F8E',
        "DMSO" : 'k', 
        "Formamide" : 'k', 
        "Water" : 'k' 
}


data = np.genfromtxt('data.txt',skip_header=1,dtype=float) 

f = open(outfile,'w') 

for i in range(1,len(data[0])) : 
    #if legend[i] == "N118X" : continue 
    plt.scatter(data[:,0],data[:,i],label=legend[i],color=nameToColorKeys[legend[i]]) 

    slope, intercept, r_value, p_value, std_error = linregress(data[:,0],data[:,i]) 
    if math.isnan(slope) : 
        slope, intercept, r_value, p_value, std_error = linregress(data[:-1,0],data[:-1,i]) 
    plt.plot(data[:,0],slope*data[:,0]+intercept,nameToColorKeys[legend[i]])


    print "%10s\t%5f\t%0.4f"%(legend[i],slope,std_error) 

    f.write("%10s\t%5f\t%0.4f\n"%(legend[i],slope,std_error) ) 

    

f.close() 

plt.legend() 
plt.savefig('ftls.png',format='png',dpi=500) 
