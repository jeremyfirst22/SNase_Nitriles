import glob 
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.colors import LogNorm
from matplotlib import rc_file
import os
from sys import exit

figCols=2
figRows=5

rcFile = 'rc_files/paper.rc'

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

rc_file(rcFile)

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

datafiles = glob.glob('B_State/*/hbond_2/persistent.xvg') 

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row',figsize=(5.3,4)) 
fig.subplots_adjust(wspace=0,hspace=0.35) 
fig.text(0.5,0.04, "Time (ns)", ha='center', va='center') 
fig.text(0.05,0.50, "Number of hydrogen bonds", rotation='vertical',ha='center', va='center') 

with open('fits/water_chosen_decay.xvg') as f: 
    lines = f.readlines() 
    fitKeys= {} 
    for line in lines : 
        splitline = line.split() 
        fitKeys[splitline[0]] = splitline[1:] 

with open('fits/protein_chosen_decay.xvg') as f: 
    lines = f.readlines() 
    protFitKeys= {} 
    for line in lines : 
        splitline = line.split() 
        protFitKeys[splitline[0]] = splitline[1:] 

xmax = 0 
for index,molec in enumerate(molecList) : 
    datafile = 'B_State/GFP_%s/hbond_2/persistent.xvg'%molec
    try : 
        data = np.genfromtxt(datafile,skip_header=24) 
        data[:,0] = data[:,0] / 1000 * 4
    except : 
        print "data failed to import for %s" %datafile
        continue 

    if (np.max(data[:,0]) > xmax) : xmax = np.max(data[:,0]) 

    
    for i in np.arange(len(data)) : 
        if data[i,5] == 0 : 
            #print "i = %i, data = %i"%(i,data[i,5]) 
            break 
    #print i 
    x = np.geomspace(data[1,0], data[i,0], 100) 
    #print data[1,0]
    #print x[0] 

    print "%s lifetimes (ps) \n\tWater: "%(datafile.split('/')[1].split('_')[1]) ,
    ax = axarr[index/figCols,index%figCols]
    if len(fitKeys[datafile.split('/')[1]] ) == 4 : 
        #print "Single exponential"
        a, A, r, adjR = fitKeys[datafile.split('/')[1]]
        a, A, r, adjR = float(a), float(A), float(r), float(adjR) 
        A = A * 1000 / 4
        y = a*np.exp(-1*A * x) 
        print "%0.1f\n"%(1/A*1000),
    elif len(fitKeys[datafile.split('/')[1]] ) == 6 : 
        #print "Double exponential"
        ax.text(0.005,1,"*",color='b',fontsize='large') 
        a, A, b, B, r, adjR = fitKeys[datafile.split('/')[1]]
        a, A, b, B, r, adjR = float(a), float(A), float(b), float(B), float(r), float(adjR) 
        A = A * 1000 / 4
        B = B * 1000 / 4
        y = a*np.exp(-1*A * x) + b* np.exp(-1* B* x) 
        print "%0.1f\t%0.1f\n"%(1/A*1000, 1/B*1000), 

    ##ax.loglog(data[:,0],data[:,4],color='tab:orange',label='Nearby water') 
    ax.loglog(data[:,0],data[:,5],color='b',linestyle='-', label='HBonding water') 
    ax.loglog(x,y,color='k',linestyle='-.') 
#    ax.semilogy(data[:,0],data[:,5],color='b',linestyle='-', label='HBonding water') 
#    ax.semilogy(x,y,color='k',linestyle='-.') 
    ax.text(.2,  500, "r = %.3f"%(r) , fontsize=10,color='b') 

    for i in np.arange(len(data)) : 
        if data[i,6] == 0 : 
            break 
    x = np.geomspace(data[1,0], data[i,0], 100) 

    print "\tProtein: " , 
    if len(protFitKeys[datafile.split('/')[1]] ) == 4 : 
        #print "Single exponential"
        a, A, r, adjR = protFitKeys[datafile.split('/')[1]]
        a, A, r, adjR = float(a), float(A), float(r), float(adjR) 
        A = A * 1000 / 4
        y = a*np.exp(-1*A * x) 
        print "%0.1f\n"%(1/A * 1000 ),
    elif len(protFitKeys[datafile.split('/')[1]] ) == 6 : 
        #print "Double exponential"
        ax.text(0.01,1,"*",color='g',fontsize='large') 
        a, A, b, B, r, adjR = protFitKeys[datafile.split('/')[1]]
        a, A, b, B, r, adjR = float(a), float(A), float(b), float(B), float(r), float(adjR) 
        A = A * 1000 / 4
        B = B * 1000 / 4
        y = a*np.exp(-1*A * x) + b* np.exp(-1* B* x) 
        print "%0.1f\t%0.1f\n"%(1/A* 1000 , 1/B* 1000), 

    #print file.split('/')[1]
    if not (datafile.split('/')[1] == "GFP_N212X" or datafile.split('/')[1] == "GFP_F114X") : 
        ax.loglog(data[:,0],data[:,6],color='g',label='Hbonding protein') 
        ax.loglog(x,y,color='k',linestyle='-.') 
#        ax.semilogy(data[:,0],data[:,6],color='g',label='Hbonding protein') 
#        ax.semilogy(x,y,color='k',linestyle='-.') 
        ax.text(.2,  150, "r = %.3f"%(r) , fontsize=10,color='g') 

    ax.set_title(datafile.split('/')[1].split('_')[1], color=nameToColorKeys[datafile.split('/')[1].split('_')[1]]) 
    #ax.set_ylim([0.8,2000]) 

for index,molec in enumerate(molecList) : 
    ax = axarr[index/figCols,index%figCols]
    ax.set_xlim(0.004  ,5) 
    ax.set_ylim(0.4,2000) 
    
fig.savefig('figures/Lifetime_hbonds.pdf',format='pdf') 
plt.close() 
