import glob 
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.colors import LogNorm
import os
import matplotlib as mpl 
from scipy.stats import linregress
from matplotlib import rc_file
from sys import exit

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

datafiles = glob.glob('SNase_*/hbond/geometry.xvg') 

left, right = 0.15, 0.95
bottom, top = 0.07, 0.95
hspace, wspace = 0.3, 0.1

fig, axarr = plt.subplots(figRows,figCols,sharex='all',sharey='all',figsize=(3.25,5)) 
fig.subplots_adjust(wspace=wspace,hspace=hspace,left=left,right=right,top=top,bottom=bottom) 
fig.text((right-left)/2+left-0.04,0.02, r"$d_{\rm{NH}}$ ($\rm{\AA}$)", ha='center', va='center') 
fig.text(0.01,(top-bottom)/2+bottom, r"$\theta_1$ (deg)", ha='left', va='center',rotation='vertical') 

xmin,xmax = 1.45, 2.45
ymin,ymax = 99,180 

for index,molec in enumerate(molecList) : 
    file="SNase_%s/hbond/geometry.xvg"%molec
    data = np.genfromtxt(file,skip_header=23) 
    try :
        data[:,0] = data[:,0] / 1000 * 4
    except IndexError : 
        print "%s is empty. Skipping"%file      
        data = np.empty((0,4),int)  
        data = np.append(data,[[0,0,0,0]],axis=0) 

    for i in range(len(data[:,3]))  : 
        if data[i,3] > 180 : 
            data[i,3] -= 180 
    ax = axarr[index%figRows,index/figRows]

    xbins, ybins = np.arange(xmin, xmax, 0.01), np.arange(ymin,ymax,1)
    z, x, y = np.histogram2d(data[:,2],data[:,3],[xbins,ybins])

    im = ax.pcolor(x,y,z.T, cmap='plasma', vmin = 0, vmax = 10) 

    ax.text(0.5,1.0,molec,color=colorDict[molec],transform=ax.transAxes,ha='center',va='bottom') 
    ax.set_ylim([ymin,179])
    ax.set_xlim([xmin,xmax]) 

fig.subplots_adjust(right=0.88) 
cbar_ax = fig.add_axes([0.90, 0.15, 0.015, 0.7]) 
fig.text(0.98,0.5, r"Counts per bin", ha='center', va='center',rotation='vertical') 
fig.colorbar(im, cax=cbar_ax) 


fig.savefig('figures/Geometries_heat.png',format='png') 
plt.close() 
#

datafiles = glob.glob('SNase_*/hbond/nw_geometry.xvg') 

fig, axarr = plt.subplots(figRows,figCols,sharex='all',sharey='all',figsize=(3.25,5)) 
fig.subplots_adjust(wspace=wspace,hspace=hspace,left=left,right=right,top=top,bottom=bottom) 
fig.text((right-left)/2+left-0.04,0.02, r"$d_{\rm{NH}}$ ($\rm{\AA}$)", ha='center', va='center') 
fig.text(0.01,(top-bottom)/2+bottom, r"$\theta_1$ (deg)", ha='left', va='center',rotation='vertical') 

xmin,xmax = 1.45, 2.45
ymin,ymax = 99,180 

for index,molec in enumerate(molecList) : 
    file="SNase_%s/hbond/nw_geometry.xvg"%molec
    data = np.genfromtxt(file,skip_header=23) 
    try :
        data[:,0] = data[:,0] / 1000 * 4
    except IndexError : 
        print "%s is empty. Skipping"%file      
        data = np.empty((0,4),int)  
        data = np.append(data,[[0,0,0,0]],axis=0) #
    for i in range(len(data[:,3]))  : 
        if data[i,3] > 180 : 
            data[i,3] -= 180 
    ax = axarr[index%figRows,index/figRows]

    xbins, ybins = np.arange(xmin, xmax, 0.01), np.arange(ymin,ymax,1)
    z, x, y = np.histogram2d(data[:,2],data[:,3],[xbins,ybins])

    im = ax.pcolor(x,y,z.T, cmap='plasma', vmin = 0, vmax = 10) 

    ax.text(0.5,1.0,molec,color=colorDict[molec],transform=ax.transAxes,ha='center',va='bottom')
    ax.set_ylim([ymin,179])
    ax.set_xlim([xmin,xmax]) 

fig.subplots_adjust(right=0.88) 
cbar_ax = fig.add_axes([0.90, 0.15, 0.015, 0.7]) 
fig.text(0.98,0.5, r"Counts per bin", ha='center', va='center',rotation='vertical') 
fig.colorbar(im, cax=cbar_ax) 

fig.savefig('figures/Geometries_protein_heat.png',format='png') 
plt.close() 


###Plot histogrammed data
#
#datafiles = glob.glob('B_State/*/fit_hbond_with_ca/dist.poly') 
#
###Theta 1 analyis
fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0.10,hspace=0.35,left=0.15,right=0.98,top=0.93,bottom=0.1) 
fig.text(0.5,0.02, r"$\theta_1$ (deg)", ha='center', va='center') 
fig.text(0.04,0.5, r"Occurences", ha='center', va='center',rotation='vertical') 

figDens, axarrDens = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
figDens.subplots_adjust(wspace=0.10,hspace=0.35,left=0.15,right=0.98,top=0.93,bottom=0.1) 
figDens.text(0.5,0.02, r"$\theta_1$ (deg)", ha='center', va='center') 
figDens.text(0.03,0.5, r"Probability density (Occurences/$\rm{\AA}^3$)", ha='center', va='center',rotation='vertical') 

f2, ax2 = plt.subplots(1,1,figsize=(4.3,3) )
f2.subplots_adjust(left=0.20,bottom=0.13,right=0.95,top=0.95) 
f2.text(0.5,0.02, r"$\theta_1$ (deg)", ha='center', va='center') 
f2.text(0.06,0.5, r"$\tilde{\nu}$ (cm$^{-1}$)", ha='center', va='center',rotation='vertical') 

r = 2.45 

absMax = {} 
with open('Exp_data/abs_data2.dat') as f : 
    for line in f : 
        if not line.startswith('#') : 
            key = line.split()[0] 
            value = float(line.split()[1]) 
            absMax[key] = value 


accumAvgTheta = [] 
accumAbsMax = [] 
for index, molec in enumerate(molecList) : 
    #print molec
    fileWhis = 'SNase_%s/fit_hbond/theta1.his'%(molec) 
    fileWpoly = 'SNase_%s/fit_hbond/theta1.poly'%(molec) 
    fileWgeom = 'SNase_%s/fit_hbond/geometry.xvg'%(molec) 
    filePhis = 'SNase_%s/fit_hbond/nw_theta1.his'%(molec) 
    filePpoly = 'SNase_%s/fit_hbond/nw_theta1.poly'%(molec) 
    filePgeom = 'SNase_%s/fit_hbond/nw_geometry.xvg'%(molec) 
    

    solvent = True 
    headlines = 0 
    with open(fileWgeom) as f :
        lines = f.readlines() 
        for line in lines : 
            if line.startswith('#') or line.startswith('@') : 
                 headlines += 1 
            else : 
                 break 
    try : 
        dataWhis = np.genfromtxt(fileWhis) 
        dataWpoly = np.genfromtxt(fileWpoly) 
        dataWgeom = np.genfromtxt(fileWgeom,skip_header=headlines) 
    except IOError : 
        dataWhis, dataWpoly = np.array([[0,0]]), np.array([[0,0]]) 
        solvent = False 

    protein = True 
    headlines = 0 
    with open(fileWgeom) as f :
        lines = f.readlines() 
        for line in lines : 
            if line.startswith('#') or line.startswith('@') : 
                 headlines += 1 
            else : 
                 break 
    try : 
        dataPhis = np.genfromtxt(filePhis) 
        dataPpoly = np.genfromtxt(filePpoly) 
        dataPgeom = np.genfromtxt(filePgeom,skip_header=headlines) 
    except IOError : 
        dataPhis, dataPpoly = np.array([[0,0]]), np.array([[0,0]]) 
        protein = False

    if solvent : 
        dataWhis[:,1] *= len(dataWgeom[:,0])  ## The histograms are unfortunately normalized. To un-normalize, 
        dataWpoly[:,1] *= len(dataWgeom[:,0])  # we mutiply each bin times the magnitude of the histogram (here),  
                                               # and the binSize (below in for loop) 
        anglesW = dataWpoly[:,0] 
        probsW = dataWpoly[:,1]

        volumesW = np.zeros_like(anglesW) 
        for i in range(len(volumesW)) : 
            if not i == len(volumesW) - 1 : 
                binSize = anglesW[i+1] - anglesW[i]
            else : 
                binSize = anglesW[i] - anglesW[i-1]
            dataWhis[i,1] *= binSize
            dataWpoly[i,1] *= binSize
            volumesW[i] = 2*np.pi * r**3 / 3 * (-np.cos((anglesW[i]+binSize/2) * np.pi / 180)  + np.cos((anglesW[i]-binSize/2)* np.pi / 180.) )
        probsW = probsW / volumesW

        avgAngle = np.average(anglesW,weights=probsW) 
        variance = np.average((anglesW- avgAngle)**2,weights=probsW) 

        avgAngleW = 0
        for i in range(len(anglesW)) : 
            avgAngleW = avgAngleW + probsW[i] * anglesW[i]

    if protein : 
        dataPhis[:,1] *= len(dataPgeom[:,0]) 
        dataPpoly[:,1] *= len(dataPgeom[:,0]) 

        anglesP = dataPpoly[:,0] 
        probsP = dataPpoly[:,1]

        volumesP = np.zeros_like(anglesP) 
        for i in range(len(volumesP)) : 
            if not i == len(volumesP) - 1 : 
                binSize = anglesP[i+1] - anglesP[i]
            else : 
                binSize = anglesP[i] - anglesP[i-1]
            dataPhis[i,1] *= binSize
            dataPpoly[i,1] *= binSize
            volumesP[i] = 2*np.pi * r**3 / 3 * (-np.cos((anglesP[i]+binSize/2) * np.pi / 180)  + np.cos((anglesP[i]-binSize/2)* np.pi / 180.) )
        #print volumesP
        probsP = probsP / volumesP

        avgAngleP = 0
        for i in range(len(anglesP)) : 
            avgAngleP = avgAngleP + probsP[i] * anglesP[i]

    if not protein and not solvent : 
        print "%s no protein or solvent hbonds. Skipping"%molec 
        continue 
    elif protein and solvent : 
        avgAngleTot = (avgAngleW + avgAngleP) / (sum(probsP)+sum(probsW))
        avgAngleP = avgAngleP / sum(probsP) 
        avgAngleW = avgAngleW / sum(probsW) 
    elif solvent : 
        avgAngleW = avgAngleW / sum(probsW) 
        avgAngleTot = avgAngleW 
    elif protein : 
        avgAngleP = avgAngleP / sum(probsP) 
        avgAngleTot = avgAngleP 
    else : 
        print "How did you get here?" 
        exit() 

    ax = axarr[index/figCols,index%figCols]
    #axD = axarrDens[index/figCols,index%figCols]

    if True  : #not molec == "A90X" and not molec == "V66X" : 
        try : 
            ax2.scatter(avgAngleTot,absMax[molec], marker='P',label=molec,color=colorDict[molec],edgecolor='none',s=100,zorder=10)  
            accumAbsMax.append(absMax[molec]) 
            accumAvgTheta.append(avgAngleTot) 
            print molec, avgAngleTot
        except KeyError : 
            print "No key found for %s"%molec
            continue 
    else : 
            ax2.scatter(avgAngleTot,absMax[molec], marker='P',label=molec,color=colorDict[molec],edgecolor='none',s=100,zorder=10,alpha=0.25)  


    #ax2.axhline(2227.5,color='k',linestyle='--') 
    #ax2.axhline(2235.9,color='#66CDFF',linestyle='--') 
    
    ax.scatter(dataWhis[:,0],dataWhis[:,1],color='b',edgecolor='none',s=2.5) 
    ax.scatter(dataPhis[:,0],dataPhis[:,1],color='g',edgecolor='none',s=2.5) 
    ax.plot(dataWpoly[:,0],dataWpoly[:,1],color='b',linestyle='--')#,linewidth=0.5) 
    ax.plot(dataPpoly[:,0],dataPpoly[:,1],color='g',linestyle='--')#,linewidth=0.5) 
    #axD.plot(anglesW,probsW,color='b',linestyle='-') 
    #axD.plot(anglesP,probsP,color='g',linestyle='-') 
#    ax.axvline(avgAngleW, color='b', linestyle='--') 
#    ax.axvline(avgAngleP, color='g', linestyle='--') 
#    ax.axvline(avgAngleTot, color='k', linestyle='--') 
    #axD.axvline(avgAngleW, color='b', linestyle='--',zorder=1) 
    #axD.axvline(avgAngleP, color='g', linestyle='--',zorder=1) 
    #axD.axvline(avgAngleTot, color='k', linestyle='--',zorder=1) 

    ax.set_title(molec,color=colorDict[molec]) 
    ax.set_ylim(0,450  ) 
    ax.set_xlim(100,180) 

    #axD.set_title(molec,color=colorDict[molec]) 
    #axD.set_ylim(0,1000 ) 
    #axD.set_xlim(100,180) 

ax2.legend(loc=2,edgecolor='k',framealpha=1.0) 
ax2.set_ylim(2155.8,2169) 
ax2.set_xlim(105,170) 

slope,intercept,r_value,p_value,std_error = linregress(accumAvgTheta, accumAbsMax)
print "r = %f, p = %f"%(r_value,p_value) 

xs = np.linspace(min(accumAvgTheta), max(accumAvgTheta) ) 
ys = slope * xs + intercept
ax2.plot(xs, ys, label="r = %0.3f"%r_value,color='k')

f2.text(0.80,0.23,r"r = %0.3f"%r_value) 

fig.savefig('figures/Geometries_angles.png',format='png') 
#figDens.savefig('figures/Geometries_angles_density.png',format='png') 
f2.savefig('figures/abs_max_vs_max_theta.pdf',format='pdf') 
plt.close() 
