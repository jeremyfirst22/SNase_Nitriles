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

#left, right = 0.15, 0.95
#bottom, top = 0.07, 0.95
#hspace, wspace = 0.3, 0.1
#
#fig, axarr = plt.subplots(figRows,figCols,sharex='all',sharey='all',figsize=(3.25,5)) 
#fig.subplots_adjust(wspace=wspace,hspace=hspace,left=left,right=right,top=top,bottom=bottom) 
#fig.text((right-left)/2+left-0.04,0.02, r"$d_{\rm{NH}}$ ($\rm{\AA}$)", ha='center', va='center') 
#fig.text(0.01,(top-bottom)/2+bottom, r"$\theta_1$ (deg)", ha='left', va='center',rotation='vertical') 
#
#xmin,xmax = 1.45, 2.45
#ymin,ymax = 99,180 
#
###########################################################################
##### Heat map
###########################################################################
#
#for index,molec in enumerate(molecList) : 
#    fileW="SNase_%s/hbond/geometry.xvg"%molec
#    fileP="SNase_%s/hbond/nw_geometry.xvg"%molec
#    dataW = np.genfromtxt(fileW,skip_header=23) 
#    dataP= np.genfromtxt(fileP,skip_header=23) 
#
#    if len(dataW) == 0 and len(dataP) == 0 :            ##no hydrogen bonds at all
#        data = np.empty((0,4),int)  
#        data = np.append(data,[[0,0,0,0]],axis=0) 
#    elif not len(dataW) == 0 and not len(dataP) == 0 :  ##both solv and prot hbonds
#        data = np.vstack((dataW,dataP)) 
#    elif len(dataW) == 0 :                              ##only prot hbonds
#        data = dataP
#    else :                                             ##only wat hbonds
#        data = dataW 
#
#    for i in range(len(data[:,3]))  : 
#        if data[i,3] > 180 : 
#            data[i,3] -= 180 
#    ax = axarr[index%figRows,index/figRows]
#
#    xbins, ybins = np.arange(xmin, xmax, 0.01), np.arange(ymin,ymax,1)
#    z, x, y = np.histogram2d(data[:,2],data[:,3],[xbins,ybins])
#
#    im = ax.pcolor(x,y,z.T, cmap='plasma', vmin = 0, vmax = 10) 
#
#    ## Overlay contour lines 
#    x = data[:,2] ; y = data[:,3]
#
#    ax.text(0.5,1.0,molec,color=colorDict[molec],transform=ax.transAxes,ha='center',va='bottom')
#    ax.set_ylim([ymin,179])
#    ax.set_xlim([xmin,xmax]) 
#
#fig.subplots_adjust(right=0.88) 
#cbar_ax = fig.add_axes([0.90, 0.15, 0.015, 0.7]) 
#fig.text(0.98,0.5, r"Counts per bin", ha='center', va='center',rotation='vertical') 
#fig.colorbar(im, cax=cbar_ax) 
#
#fig.savefig('figures/Geometries_combined_heat.png',format='png') 
#plt.close() 
##
#
###########################################################################
##### Theta 1 
###########################################################################
#fig, axarr = plt.subplots(figRows,figCols,sharex='all',sharey='all',figsize=(3,5) )
#fig.subplots_adjust(wspace=0.15,hspace=0.35,left=0.15,right=0.95,top=0.93,bottom=0.1) 
#fig.text(0.5,0.02, r"Time", ha='center', va='center') 
#fig.text(0.01,0.5, r"$\theta_1$ (deg)", ha='left', va='center',rotation='vertical') 
#
#fig2, ax2 = plt.subplots(1,1,sharex='col',sharey='row') 
#fig2.subplots_adjust(wspace=0.10,hspace=0.35,left=0.12,right=0.98,top=0.93,bottom=0.1) 
#fig2.text(0.5,0.02, r"$\theta_1$ (deg)", ha='center', va='center') 
#fig2.text(0.01,0.5, r"Occurances (arb. units)", ha='left', va='center',rotation='vertical') 
#
#for index,molec in enumerate(molecList) : 
#    fileW="SNase_%s/hbond/geometry.xvg"%molec
#    fileP="SNase_%s/hbond/nw_geometry.xvg"%molec
#    dataW = np.genfromtxt(fileW,skip_header=23) 
#    dataP= np.genfromtxt(fileP,skip_header=23) 
#
#    ax = axarr[index%figRows,index/figRows]
#
#    if not len(dataW) == 0 : 
#        for i in range(len(dataW[:,3])) : 
#            if dataW[i,3] > 180 : 
#                dataW[i,3] -= 180 
#        dataW[:,0] = dataW[:,0] * 4 / 1000. 
#        x = dataW[:,0] ; y = dataW[:,3]
#        ax.scatter(x,y,s=0.1,color='b') 
#    if not len(dataP) == 0 : 
#        for i in range(len(dataP[:,3])) : 
#            if dataP[i,3] > 180 : 
#                dataP[i,3] -= 180 
#        dataP[:,0] = dataP[:,0] * 4 / 1000. 
#        x = dataP[:,0] ; y = dataP[:,3]
#        ax.scatter(x,y,s=0.1,color='g') 
#
#    ax.set_title(molec,color=colorDict[molec]) 
#    ax.set_ylim([100,179])
#    ax.set_xlim([0,100]) 
#
#    fileHis = "SNase_%s/fit_hbond/theta1.his"%molec
#    filePoly = "SNase_%s/fit_hbond/theta1.poly"%molec
#    
#    dataHis = np.genfromtxt(fileHis) 
#    dataPoly= np.genfromtxt(filePoly) 
#
#    dataHis[:,1] *= (len(dataW) + len(dataP)) / 1000 
#    dataPoly[:,1] *= (len(dataW) + len(dataP)) / 1000 
#
#    ax2.scatter(dataHis[:,0],dataHis[:,1],s=0.1,color=colorDict[molec])
#    ax2.plot(dataPoly[:,0],dataPoly[:,1],color=colorDict[molec]) 
#
#fig.savefig('figures/Geometries_time_vs_theta1.png',format='png') 
#fig2.savefig('figures/Geometries_hist_theta1.png',format='png') 
#plt.close() 
#
###########################################################################
##### Theta 2 
###########################################################################
#fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
#fig.subplots_adjust(wspace=0.10,hspace=0.35,left=0.12,right=0.98,top=0.93,bottom=0.1) 
#fig.text(0.5,0.02, r"Time", ha='center', va='center') 
#fig.text(0.01,0.5, r"$\theta_2$ (deg)", ha='left', va='center',rotation='vertical') 
#
#fig2, ax2 = plt.subplots(1,1,sharex='col',sharey='row') 
#fig2.subplots_adjust(wspace=0.10,hspace=0.35,left=0.12,right=0.98,top=0.93,bottom=0.1) 
#fig2.text(0.5,0.02, r"$\theta_2$ (deg)", ha='center', va='center') 
#fig2.text(0.01,0.5, r"Occurances (arb. units)", ha='left', va='center',rotation='vertical') 
#
#for index,molec in enumerate(molecList) : 
#    fileW="SNase_%s/hbond/geometry.xvg"%molec
#    fileP="SNase_%s/hbond/nw_geometry.xvg"%molec
#    dataW = np.genfromtxt(fileW,skip_header=23) 
#    dataP= np.genfromtxt(fileP,skip_header=23) 
#
#    ax = axarr[index/figCols,index%figCols]
#
#    if not len(dataW) == 0 : 
#        for i in range(len(dataW[:,3])) : 
#            if dataW[i,3] > 180 : 
#                dataW[i,3] -= 180 
#        dataW[:,0] = dataW[:,0] * 4 / 1000. 
#        x = dataW[:,0] ; y = dataW[:,4]
#        ax.scatter(x,y,s=0.1,color='b') 
#    if not len(dataP) == 0 : 
#        for i in range(len(dataP[:,3])) : 
#            if dataP[i,3] > 180 : 
#                dataP[i,3] -= 180 
#        dataP[:,0] = dataP[:,0] * 4 / 1000. 
#        x = dataP[:,0] ; y = dataP[:,4]
#        ax.scatter(x,y,s=0.1,color='g') 
#
#    ax.set_title(molec,color=colorDict[molec]) 
#    ax.set_ylim([100,179])
#    ax.set_xlim([0,100]) 
#
#    fileHis = "SNase_%s/fit_hbond/theta2.his"%molec
#    filePoly = "SNase_%s/fit_hbond/theta2.poly"%molec
#    
#    dataHis = np.genfromtxt(fileHis) 
#    dataPoly= np.genfromtxt(filePoly) 
#
#    dataHis[:,1] *= (len(dataW) + len(dataP)) / 1000 
#    dataPoly[:,1] *= (len(dataW) + len(dataP)) / 1000 
#
#    ax2.scatter(dataHis[:,0],dataHis[:,1],s=0.1,color=colorDict[molec])
#    ax2.plot(dataPoly[:,0],dataPoly[:,1],color=colorDict[molec]) 
#
#fig.savefig('figures/Geometries_time_vs_theta2.png',format='png') 
#fig2.savefig('figures/Geometries_hist_theta2.png',format='png') 
#plt.close() 
#
###########################################################################
##### Distance N-H
###########################################################################
#fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
#fig.subplots_adjust(wspace=0.10,hspace=0.35,left=0.12,right=0.98,top=0.93,bottom=0.1) 
#fig.text(0.5,0.02, r"Time", ha='center', va='center') 
#fig.text(0.01,0.5, r"$d_{\rm{NH}}$ ($\AA$)", ha='left', va='center',rotation='vertical') 
#
#fig2, ax2 = plt.subplots(1,1,sharex='col',sharey='row') 
#fig2.subplots_adjust(wspace=0.10,hspace=0.35,left=0.12,right=0.98,top=0.93,bottom=0.1) 
#fig2.text(0.5,0.02, r"$d_{\rm{NH}}$ ($\AA$)", ha='center', va='center') 
#fig2.text(0.01,0.5, r"Occurances (arb. units)", ha='left', va='center',rotation='vertical') 
#
#for index,molec in enumerate(molecList) : 
#    fileW="SNase_%s/hbond/geometry.xvg"%molec
#    fileP="SNase_%s/hbond/nw_geometry.xvg"%molec
#    dataW = np.genfromtxt(fileW,skip_header=23) 
#    dataP= np.genfromtxt(fileP,skip_header=23) 
#
#    ax = axarr[index/figCols,index%figCols]
#
#    if not len(dataW) == 0 : 
#        dataW[:,0] = dataW[:,0] * 4 / 1000. 
#        x = dataW[:,0] ; y = dataW[:,2]
#        ax.scatter(x,y,s=0.1,color='b') 
#    if not len(dataP) == 0 : 
#        dataP[:,0] = dataP[:,0] * 4 / 1000. 
#        x = dataP[:,0] ; y = dataP[:,2]
#        ax.scatter(x,y,s=0.1,color='g') 
#
#    ax.set_title(molec,color=colorDict[molec]) 
#    ax.set_ylim([1.5,2.5])
#    ax.set_xlim([0,100]) 
#
#    fileHis = "SNase_%s/fit_hbond/dist.his"%molec
#    filePoly = "SNase_%s/fit_hbond/dist.poly"%molec
#    
#    dataHis = np.genfromtxt(fileHis) 
#    dataPoly= np.genfromtxt(filePoly) 
#
#    dataHis[:,1] *= (len(dataW) + len(dataP)) / 1000 
#    dataPoly[:,1] *= (len(dataW) + len(dataP)) / 1000 
#
#    ax2.scatter(dataHis[:,0],dataHis[:,1],s=0.1,color=colorDict[molec])
#    ax2.plot(dataPoly[:,0],dataPoly[:,1],color=colorDict[molec]) 
#
#    print "%s\tWater: %i\tProtein: %i"%(molec,len(dataW),len(dataP)) 
#
#fig.savefig('figures/Geometries_time_vs_dist.png',format='png') 
#fig2.savefig('figures/Geometries_hist_dist.png',format='png') 
#plt.close() 

##########################################################################
#### Lifetime of hbonds 
##########################################################################

def exp_decay(x, a,b) : 
    return a*np.exp(-b*x)
def double_exp_decay(x, a,b, c,d) : 
    return a*np.exp(-b*x) + c*np.exp(-d*x)
def three_exp_decay(x, a,b, c,d, e,f) : 
    return a*np.exp(-b*x) + c*np.exp(-d*x) + e*np.exp(-f*x)

left, right = 0.15, 0.95
bottom, top = 0.08, 0.92
wspace,hspace = 0.1,0.3

fig, axarr = plt.subplots(figRows,figCols,sharex='all',sharey='all',figsize=(3.25,6)) 
fig.subplots_adjust(wspace=wspace,hspace=hspace,bottom=bottom,top=top,right=right,left=left)
fig.text(left+(right-left)/2,0.01, r"$t$ (ns)", ha='center', va='bottom') 
fig.text(0.01,bottom+(top-bottom)/2, r"Number of hydrogen bonds that persist for a least time $t$", ha='left', va='center',rotation='vertical') 

for index,molec in enumerate(molecList) : 

    ax = axarr[index%figRows,index/figRows]

    file ="SNase_%s/hbond_2/persistent.xvg"%molec
    data = np.genfromtxt(file,skip_header=24) 

    x = data[:,0] 
    y = data[:,5]

    #x /= 1000 ##ps -> ns 

    amin = 0 
    amax = np.max(y)
    bmin = 0 
    bmax = 40 


    if molec in ["A109X","V66X","A58X","V23X"] : 
        numPeaks, func = 2, double_exp_decay
        ax.text(0.98,0.98,"**",transform=ax.transAxes,va='top',ha='right')
    elif molec in ["N118X"] :
        numPeaks, func = 3, three_exp_decay
        ax.text(0.98,0.98,"***",transform=ax.transAxes,va='top',ha='right')
    else : 
        numPeaks, func = 1, exp_decay
        ax.text(0.98,0.98,"*",transform=ax.transAxes,va='top',ha='right')

    mins = np.tile([amin,bmin],numPeaks)
    maxs = np.tile([amax,bmax],numPeaks)
    bounds = (mins,maxs) 
    popt, pcov = curve_fit(func, x, y,bounds=bounds,maxfev=100000) 

    xs = np.linspace(np.min(x),np.max(x),len(x)*10) 
    fit = func(xs,*popt)
    ax.plot(xs,fit,'k--',alpha=0.75) 

    #print molec, popt
    print "%10s"%molec, 
    for item in popt[1::2] : 
        print "\t%8.3f"%item,
    print

    #for i in range(len(popt)/2) : 
    #    params = popt[i*2:i*2+2]
    #    fit = exp_decay(xs, *params) 
    #    ax.plot(xs,fit,'k--',alpha=0.5) 


    ax.plot(x,y,color='b')
    #ax.plot(data[:,0],data[:,6],color='g')


    ax.set_title(molec,color=colorDict[molec]) 
    ax.set_ylim([1,1000])
    ax.set_xlim([1,100]) 
    #ax.set_yscale('log') 
    #ax.set_xscale('log') 

fig.savefig('figures/Lifetimes.png',format='png') 
plt.close() 
