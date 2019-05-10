import numpy as np 
import matplotlib.pyplot as plt 
import math
from scipy.stats import linregress 
from matplotlib import rc_file
from sys import exit

#sys.path.insert(0, "/Users/jfirst/
from adjustText import adjust_text

file = 'Exp_data/FTLS.dat'
outfile = 'Exp_data/ftls_fits.dat'

tempList = [] 
ftlsDict = {} 
with open(file) as f : 
    for item in f.readline().split() : 
        if "std" in item : continue 
        if "#" in item : continue 
        if "deg" in item : continue 
        else : 
            tempList.append(float(item)) 

    for line in f.readlines() : 
        molec = line.split()[0]
        peaks = line.split()[1::2]
        stds  = line.split()[2::2]
        peaks = np.array(peaks, dtype=float) 
        stds = np.array(stds, dtype=float) 
        ftlsDict[molec] = [peaks,stds]

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
    "A58X",
    "Water",
    "DMSO"
]

rc_file('rc_files/paper.rc') 

left, right = 0.17, 0.98
bottom, top = 0.15, 0.98

fig, ax    = plt.subplots(1,1,figsize=(3.42,2.25)) 
fig.subplots_adjust(left=left, bottom=bottom,right=right,top=top) 
fig.text((right-left)/2+left,0.01,           r"Temperature ($^{\circ}$C)",ha='center', va='bottom')
fig.text(0.01,(top-bottom)/2+bottom,         r"Nitrile frequency shift (cm$^{-1}$)                ",ha='left',va='center',rotation='vertical')

f = open(outfile,'w') 

texts = [] 
#for i in range(1,len(data[0])) : 
for molec in molecList : 
    temps = np.array(tempList,dtype=float) 
    try : 
        peaks = ftlsDict[molec][0]
        stds = ftlsDict[molec][1]
    except KeyError : 
        print "No data found for key %s"%molec
        continue 
    assert len(peaks) == len(temps) == len(stds) 

    if molec == "DMSO" : 
        peaks = peaks[1:]
        stds = stds[1:]
        temps = temps[1:]

        slope, intercept, r_value, p_value, std_error = linregress(temps,peaks)
        peaks -= slope*tempList[0] + intercept 
    else: 
        peaks -= peaks[0]

    for i,temp in enumerate(temps) : 
        ax.scatter(temp,peaks[i],color=colorDict[molec]) 
        #ax.errorbar(temp,peaks[i],yerr=stds[i],color=colorDict[molec]) 
    #ax.scatter(data[:,0],data[:,i],label=legend[i],color=colorDict[legend[i]]) 

    #if not molec in ["NX118X"] : 
    #    slope, intercept, r_value, p_value, std_error = linregress(temps,peaks) 
    #else : 
    #    slope, intercept, r_value, p_value, std_error = linregress(np.append(temps[0],temps[2:]),np.append(peaks[0],peaks[2:])) 
    #if math.isnan(slope) : 
    #    slope, intercept, r_value, p_value, std_error = linregress(data[:-1,0],data[:-1,i]) 
    slope, intercept, r_value, p_value, std_error = linregress(temps,peaks) 
    lines = ax.plot(temps,slope*temps+intercept,colorDict[molec]) 

    x = temps[-1]
    y = slope*x+intercept
    x += 6.30
    if molec == "V66X" : 
        y += 0.22
    elif molec == "V23X" : 
        y += 0.15
    elif molec == "L25X" : 
        y += 0.07
    elif molec == "DMSO" : 
        y -= 0.01 
    elif molec == "A90X" : 
        y -= 0.04
    elif molec == "T62X" : 
        y -= 0.06
    if molec == "I92X" : 
        y -= 0.09 
    if molec == "A109X" : 
        y -= 0.09 
    elif molec == "A58X" : 
        y -= 0.10
    texts.append(ax.text(x,y,molec,color=colorDict[molec], ha='right',va='center') ) #,color=colorDict[molec],va='center',ha='left') )


    print "%10s\t%5f\t%0.4f"%(molec,slope,std_error) 

    f.write("%10s\t%5f\t%0.4f\n"%(molec,slope,std_error) ) 

f.close() 

ax.set_xlim(4,42)

#adjust_text(texts\
#        #,add_objects=lines
#        ,autoalign='y'\
#        ,expand_text=(0.0,1.0)\
#        ,only_move={'points':'y', 'text':'xy', 'objects':'y'}\
#        #,save_steps=True,save_prefix="step"\
#        )

#ax.legend(loc=3) 
fig.savefig('figures/ftls.png',format='png',dpi=500) 
