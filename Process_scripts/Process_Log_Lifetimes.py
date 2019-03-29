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

def exp_decay(x, a,k) : 
    return np.log(a*np.exp(-k*x))
    #return a*x + b 
def double_exp_decay(x, a1,k1, a2,k2) : 
    return np.log(a1*np.exp(-k1*x) + a2*np.exp(-k2*x))
def three_exp_decay(x, a1,k1, a2,k2, a3,k3) : 
    return np.log(a1*np.exp(-k1*x) + a2*np.exp(-k2*x) + a3*np.exp(-k3*x))
#def double_exp_decay(x, a,b, c,d) : 
    #return a*np.exp(-1/b*x) + c*np.exp(-1/d*x)
#    return (-1/b*x + a )*(-1/d*x + c)

left, right = 0.15, 0.95
bottom, top = 0.06, 0.98
wspace,hspace = 0.1,0.3

fig, axarr = plt.subplots(figRows,figCols,sharex='none',sharey='all',figsize=(3.25,6)) 
fig.subplots_adjust(wspace=wspace,hspace=hspace,bottom=bottom,top=top,right=right,left=left)
fig.text(left+(right-left)/2,0.01, r"$t$ (ns)", ha='center', va='bottom') 
fig.text(0.01,bottom+(top-bottom)/2, r"Number of hydrogen bonds that persist for a least time $t$", ha='left', va='center',rotation='vertical') 

fig2, axarr2 = plt.subplots(figRows,figCols,sharex='none',sharey='all',figsize=(3.25,6)) 
fig2.subplots_adjust(wspace=wspace,hspace=hspace,bottom=bottom,top=top,right=right,left=left)
fig2.text(left+(right-left)/2,0.01, r"$t$ (ns)", ha='center', va='bottom') 
fig2.text(0.01,bottom+(top-bottom)/2, r"Number of hydrogen bonds that persist for a least time $t$", ha='left', va='center',rotation='vertical') 

print "%10s   %8s  %8s   %8s  %8s   %8s  %8s"%("molec","a1","k1","a2","k2","a3","k3")
for index,molec in enumerate(molecList) : 
    ax = axarr[index%figRows,index/figRows]
    ax2= axarr2[index%figRows,index/figRows]


    ax.text(0.10,0.90,molec,color=colorDict[molec],ha='left',va='top',transform=ax.transAxes)
    ax.set_ylim([np.log(1),np.log(5000)])
    ax.set_yticks(np.arange(np.log(1),np.log(5000),2) )
    ax.set_yticklabels([1,"10$^2$","10$^4$","10$^6$","10$^8$"]) 
    ax.set_xlim([0.004,3.25]) 
    ax.set_xscale('log') 

    ax2.text(0.10,0.90,molec,color=colorDict[molec],ha='left',va='top',transform=ax.transAxes)
    ax2.set_ylim([np.log(1),np.log(5000)])
    ax2.set_yticks(np.arange(np.log(1),np.log(5000),2) )
    ax2.set_yticklabels([1,"10$^2$","10$^4$","10$^6$","10$^8$"]) 
    ax2.set_xlim([0.004,3.25]) 

    file ="SNase_%s/hbond_2/persistent.xvg"%molec
    data = np.genfromtxt(file,skip_header=24) 

    x = data[:,0] 
    x *= 4 ##frames -> ps 
    x /= 1000 ##ps -> ns 

    y = data[:,5]
    y = np.trim_zeros(y)

    x = x[:len(y)]
    ax.scatter(x,np.log(y),color='b',s=2)
    ax2.scatter(x,np.log(y),color='b',s=2)

    #y = np.trim_zeros(y-1)+1 
    #y = np.trim_zeros(y-2)+2 
    #y = np.trim_zeros(y-3)+3 
    #y = np.trim_zeros(y-4)+4 

    x = x[:len(y)]
    y = np.log(y) 




    if False :#:molec in ["A90X","V66X","A109X","V23X","L25X","A58X"] : 
        numPeaks, func = 1, exp_decay
        #ax.text(0.98,0.98,"*",transform=ax.transAxes,va='top',ha='right')
        #ax2.text(0.98,0.98,"*",transform=ax.transAxes,va='top',ha='right')
    elif False : #molec in ["T62X","I92X","N118X"] : 
        numPeaks, func = 2, double_exp_decay
        #ax.text(0.98,0.98,"**",transform=ax.transAxes,va='top',ha='right')
        #ax2.text(0.98,0.98,"**",transform=ax.transAxes,va='top',ha='right')
    elif True : #molec in ["A109X","N118X"] : 
        numPeaks, func = 3, three_exp_decay
        #ax.text(0.98,0.98,"***",transform=ax.transAxes,va='top',ha='right')
        #ax2.text(0.98,0.98,"***",transform=ax.transAxes,va='top',ha='right')
    else : 
        #print "How did you get here? %s"%molec
        #ax.text(0.98,0.85,"No fit",transform=ax.transAxes,va='top',ha='right')
        #ax2.text(0.98,0.85,"No fit",transform=ax.transAxes,va='top',ha='right')
        continue 

    amin = 0.000 
    amax = np.max(np.exp(y))*10
    kmin =  0.0004
    kmax =  250

    mins = np.tile([amin,kmin],numPeaks)
    maxs = np.tile([amax,kmax],numPeaks)
    bounds = (mins,maxs) 
    
    #print molec,"\t",bounds

    popt, pcov = curve_fit(func, x, y,bounds=bounds,maxfev=5000) 

    xs = np.linspace(np.min(x)+0.004,np.max(x),len(x)*10) 
    fit = func(xs,*popt)
    ax.plot(xs,fit,'k--',zorder=3) 

    xs = np.linspace(0,3,100)
    fit = func(xs,*popt)
    ax2.plot(xs,fit,'k--',zorder=3) 

    fit = func(x,*popt)
    residuals = y - fit
    ss_res = np.sum(residuals**2) 
    ss_tot = np.sum((y-np.mean(y))**2) 

    r_squared = 1 - (ss_res / ss_tot)

    ax.text(0.98,0.85,"r = %.3f"%(np.sqrt(r_squared)),transform=ax.transAxes,va='top',ha='right')
    ax2.text(0.98,0.85,"r = %.3f"%(np.sqrt(r_squared)),transform=ax.transAxes,va='top',ha='right')

    print "%10s"%molec, 
    #for item in popt[1::2] : 
    #    print "   %8.0f"%(1/item*1000),
    #print "\tR^2 = %0.3f"%((r_squared))
    for i in range(0,len(popt),2) : 
        print "   %8.1f %8.3f"%(popt[i],1/popt[i+1]*1000),
    print
    #for i in range(len(popt)/2) : 
    #    params = popt[i*2:i*2+2]
    #    fit = exp_decay(xs, *params) 
    #    ax.plot(xs,fit,'k--',alpha=0.5) 


    #ax.plot(data[:,0],data[:,6],color='g')


    #ax.set_title(molec,color=colorDict[molec]) 

fig.savefig('figures/Log_Lifetimes.png',format='png') 
fig2.savefig('figures/Lifetimes.png',format='png') 
plt.close() 
