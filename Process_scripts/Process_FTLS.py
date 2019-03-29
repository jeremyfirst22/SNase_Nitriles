import numpy as np 
import matplotlib.pyplot as plt 
import math
from scipy.stats import linregress 
from matplotlib import rc_file

#sys.path.insert(0, "/Users/jfirst/
from adjustText import adjust_text

file = 'Exp_data/FTLS/data.txt'
outfile = 'Exp_data/ftls.dat'

with open(file) as f : 
    legend = f.readlines()[0]
    legend = legend.split() 

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

data = np.genfromtxt(file,skip_header=1,dtype=float) 

rc_file('rc_files/paper.rc') 


left, right = 0.18, 0.98
bottom, top = 0.10, 0.98

fig, ax    = plt.subplots(1,1,figsize=(3.25,3.25)) 
fig.subplots_adjust(left=left, bottom=bottom,right=right,top=top) 
fig.text((right-left)/2+left,0.01,           r"Temperature ($^{\circ}$C)",ha='center', va='bottom')
fig.text(0.01,(top-bottom)/2+bottom,         r"Nitrile frequency shift (cm$^{-1}$)                ",ha='left',va='center',rotation='vertical')

f = open(outfile,'w') 

texts = [] 
for i in range(1,len(data[0])) : 
#    if legend[i] == "N118X" : continue 
    if legend[i] == "Formamide" : continue 
    ax.scatter(data[:,0],data[:,i],label=legend[i],color=colorDict[legend[i]]) 

    slope, intercept, r_value, p_value, std_error = linregress(data[:,0],data[:,i]) 
    if math.isnan(slope) : 
        slope, intercept, r_value, p_value, std_error = linregress(data[:-1,0],data[:-1,i]) 
    lines = ax.plot(data[:,0],slope*data[:,0]+intercept,colorDict[legend[i]]) 

    x = data[-1,0] 
    y = slope*x+intercept
    x += 6.5
    if legend[i] == "I92X" : 
        y -= 0.06 
    elif legend[i] == "A90X" : 
        y -= 0.07
    elif legend[i] == "V66X" : 
        y -= 0.02
    texts.append(ax.text(x,y,legend[i],color=colorDict[legend[i]], ha='right',va='center') ) #,color=colorDict[legend[i]],va='center',ha='left') )


    print "%10s\t%5f\t%0.4f"%(legend[i],slope,std_error) 

    f.write("%10s\t%5f\t%0.4f\n"%(legend[i],slope,std_error) ) 

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
