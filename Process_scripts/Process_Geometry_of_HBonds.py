import glob 
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.colors import LogNorm
import os

figCols=2
figRows=3

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

datafiles = glob.glob('*/hbond/geometry.xvg') 

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0.1) 
fig.subplots_adjust(hspace=0.25) 
fig.text(0.5,0.04, "Time (ns)", ha='center', va='center') 
fig.text(0.05,0.5, "HBond Angle", ha='center', va='center',rotation='vertical') 

for index,file in enumerate(datafiles) : 
    molec = file.split('/')[0]
    title=molec.split('_')[1]
    if title[-1] == 'X' : 
        title = title[:-1]
        title+='C$_{\\rm{SCN}}$'
 
    data = np.genfromtxt(file,skip_header=23) 
    data[:,0] = data[:,0] / 1000 * 4

    for i in range(len(data[:,2]))  : 
        if data[i,3] > 180 : 
            data[i,3] -= 180 
    ax = axarr[index/figCols,index%figCols]

    ax.scatter(data[:,0],data[:,3],s=0.1) 
    ax.set_title(title) 
    ax.set_xlim([0,50])
    ax.set_ylim([90,190])

fig.savefig('figures/Geometries_angles.png',format='png') 
plt.close() 


fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0.1) 
fig.subplots_adjust(hspace=0.25) 
fig.text(0.5,0.04, "Time (ns)", ha='center', va='center') 
fig.text(0.05,0.5, r"HBond Length ($\AA$)", ha='center', va='center',rotation='vertical') 

for index,file in enumerate(datafiles) : 
    molec = file.split('/')[0]
    title=molec.split('_')[1]
    if title[-1] == 'X' : 
        title = title[:-1]
        title+='C$_{\\rm{SCN}}$'

    data = np.genfromtxt(file,skip_header=23) 
    data[:,0] = data[:,0] / 1000 * 4

    for i in range(len(data[:,3]))  : 
        if data[i,3] > 180 : 
            data[i,3] -= 180 
    ax = axarr[index/figCols,index%figCols]

    ax.scatter(data[:,0],data[:,2],s=0.1) 
    ax.set_title(title) 
    ax.set_xlim([0,50])
    ax.set_ylim([1.4,2.5])

fig.savefig('figures/Geometries_length.png',format='png') 
plt.close()

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0.1) 
fig.subplots_adjust(hspace=0.25) 
fig.text(0.5,0.04, r"HBond Length ($\AA$)", ha='center', va='center') 
fig.text(0.05,0.5, r"HBond Angle (deg)", ha='center', va='center',rotation='vertical') 

for index,file in enumerate(datafiles) : 
    molec = file.split('/')[0]
    title=molec.split('_')[1]
    if title[-1] == 'X' : 
        title = title[:-1]
        title+='C$_{\\rm{SCN}}$'

    data = np.genfromtxt(file,skip_header=23) 
    data[:,0] = data[:,0] / 1000 * 4

    for i in range(len(data[:,3]))  : 
        if data[i,3] > 180 : 
            data[i,3] -= 180 
    ax = axarr[index/figCols,index%figCols]

    print file.split('/')[1], "\t", len(data[:,3]) 
    ax.scatter(data[:,2],data[:,3],s=0.05) 
    ax.set_title(title) 
    ax.set_ylim([90,190])
    ax.set_xlim([1.4,2.5])
    print file.split('/')[1],"\t",len(data[:,3]) 

fig.savefig('figures/Geometries_2D.png',format='png') 
plt.close() 

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0.1) 
fig.subplots_adjust(hspace=0.25) 
fig.text(0.5,0.04, r"HBond Length ($\AA$)", ha='center', va='center') 
fig.text(0.05,0.5, r"HBond Angle (deg)", ha='center', va='center',rotation='vertical') 

for index,file in enumerate(datafiles) : 
    molec = file.split('/')[0]
    title=molec.split('_')[1]
    if title[-1] == 'X' : 
        title = title[:-1]
        title+='C$_{\\rm{SCN}}$'

    data = np.genfromtxt(file,skip_header=23) 
    data[:,0] = data[:,0] / 1000 * 4

    for i in range(len(data[:,3]))  : 
        if data[i,3] > 180 : 
            data[i,3] -= 180 
    ax = axarr[index/figCols,index%figCols]

    x = data[:,2] ; y = data[:,3]
    counts,  xbins,  ybins  = np.histogram2d(x, y) #,bins=(64,64)) 

    ax.contour(counts.transpose(),extent=[xbins.min(),xbins.max(),
       ybins.min(),ybins.max()],linewidths=2,colors='black',
           linestyles='solid')

    ax.set_title(title) 
    ax.set_ylim([90,190])
    ax.set_xlim([1.4,2.5])

fig.savefig('figures/Geometries_contour.png',format='png') 
plt.close() 

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0.1) 
fig.subplots_adjust(hspace=0.25) 
fig.text(0.5,0.04, r"HBond Length ($\AA$)", ha='center', va='center') 
fig.text(0.05,0.5, r"HBond Angle (deg)", ha='center', va='center',rotation='vertical') 

for index,file in enumerate(datafiles) : 
    molec = file.split('/')[0]
    title=molec.split('_')[1]
    if title[-1] == 'X' : 
        title = title[:-1]
        title+='C$_{\\rm{SCN}}$'

    data = np.genfromtxt(file,skip_header=23) 
    data[:,0] = data[:,0] / 1000 * 4

    for i in range(len(data[:,3]))  : 
        if data[i,3] > 180 : 
            data[i,3] -= 180 
    ax = axarr[index/figCols,index%figCols]

    xbins, ybins = np.arange(1.45,2.5,0.05), np.arange(99,180,5)
    z, x, y = np.histogram2d(data[:,2],data[:,3],[xbins,ybins])

    ax.pcolor(x,y,z.T, vmin = 0, vmax = 150) 

    ## Overlay contour lines 
    x = data[:,2] ; y = data[:,3]
    counts,  xbins,  ybins  = np.histogram2d(x, y) #,bins=(64,64)) 

    ax.contour(counts.transpose(),extent=[xbins.min(),xbins.max(),
       ybins.min(),ybins.max()],linewidths=1,colors='black',
           linestyles='solid')

    ax.set_title(title) 
    ax.set_ylim([90,190])
    ax.set_xlim([1.4,2.5])


fig.savefig('figures/Geometries_heat.png',format='png') 
plt.close() 
