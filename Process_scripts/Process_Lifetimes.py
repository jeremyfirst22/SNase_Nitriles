import glob 
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.colors import LogNorm
import os
import matplotlib as mpl 
from scipy.stats import linregress
from matplotlib import rc_file
from sys import exit
from scipy.optimize import curve_fit

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

figCols=2
figRows=5

rc_file(rcFile) 

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

##########################################################################
#### Lifetime of hbonds 
##########################################################################

def exp_decay(x, a,b) : 
    return a*np.exp(-1/b*x)
def double_exp_decay(x, a,b, c,d) : 
    return a*np.exp(-1/b*x) + c*np.exp(-1/d*x)
def three_exp_decay(x, a,b, c,d, e,f) : 
    return a*np.exp(-1/b*x) + c*np.exp(-1/d*x) + e*np.exp(-1/f*x)
def four_exp_decay(x, a,b, c,d, e,f, g,h) : 
    return a*np.exp(-1/b*x) + c*np.exp(-1/d*x) + e*np.exp(-1/f*x) + g*np.exp(-1/h*x)
def five_exp_decay(x, a,b, c,d, e,f, g,h, i,j) : 
    return a*np.exp(-1/b*x) + c*np.exp(-1/d*x) + e*np.exp(-1/f*x) + g*np.exp(-1/h*x) + i*np.exp(-j*x)

left, right = 0.15, 0.95
bottom, top = 0.06, 0.98
wspace,hspace = 0.1,0.3

fig, axarr = plt.subplots(figRows,figCols,sharex='none',sharey='all',figsize=(3.25,6)) 
fig.subplots_adjust(wspace=wspace,hspace=hspace,bottom=bottom,top=top,right=right,left=left)
fig.text(left+(right-left)/2,0.01, r"$t$ (ns)", ha='center', va='bottom') 
fig.text(0.01,bottom+(top-bottom)/2, r"Number of hydrogen bonds that persist for a least time $t$", ha='left', va='center',rotation='vertical') 

for index,molec in enumerate(molecList) : 

    ax = axarr[index%figRows,index/figRows]

    file ="SNase_%s/hbond_2/persistent.xvg"%molec
    data = np.genfromtxt(file,skip_header=24) 

    x = data[:,0] 
    y = data[:,5]

    y = np.trim_zeros(y)
    x = x[:len(y)]

    x *= 4 ##frames -> ps 
    x /= 1000 ##ps -> ns 

    amin = 0 
    amax = np.max(y)
    bmin = 0.001
    bmax = 10 


    if molec in ["V66X","T62X","L38X","A58X"] : 
        numPeaks, func = 1, exp_decay
        ax.text(0.98,0.98,"*",transform=ax.transAxes,va='top',ha='right')
    elif molec in ["A90X","V23X","L25X","I92X"] : 
        numPeaks, func = 2, double_exp_decay
        ax.text(0.98,0.98,"**",transform=ax.transAxes,va='top',ha='right')
    elif molec in ["A109X","N118X"] :
        numPeaks, func = 3, three_exp_decay
        ax.text(0.98,0.98,"***",transform=ax.transAxes,va='top',ha='right')
    else : 
        print "How did you get here? %s"%molec
        continue 
    #else : 
    #    numPeaks, func = 1, exp_decay
    #    ax.text(0.98,0.98,"*",transform=ax.transAxes,va='top',ha='right')

    mins = np.tile([amin,bmin],numPeaks)
    maxs = np.tile([amax,bmax],numPeaks)
    bounds = (mins,maxs) 
    popt, pcov = curve_fit(func, x, y,bounds=bounds,maxfev=10000) 

    xs = np.linspace(np.min(x)+0.004,np.max(x),len(x)*10) 
    fit = func(xs,*popt)
    ax.plot(xs,fit,'k--') 

    #print molec, popt
    print "%10s"%molec, 
    for item in popt[1::2] : 
        print "\t%8.3f"%(item*1000),
    print

    #for i in range(len(popt)/2) : 
    #    params = popt[i*2:i*2+2]
    #    fit = exp_decay(xs, *params) 
    #    ax.plot(xs,fit,'k--',alpha=0.5) 


    ax.plot(x,y,color='b')
    #ax.plot(data[:,0],data[:,6],color='g')


    #ax.set_title(molec,color=colorDict[molec]) 
    ax.text(0.10,0.90,molec,color=colorDict[molec],ha='left',va='top',transform=ax.transAxes)
    ax.set_ylim([1,5000])
    ax.set_xlim([0.004,1]) 
    ax.set_yscale('log') 
    ax.set_xscale('log') 

fig.savefig('figures/Lifetimes.png',format='png') 
plt.close() 
