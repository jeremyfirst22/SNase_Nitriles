import glob 
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.colors import LogNorm
import os
import matplotlib as mpl 
from scipy.stats import linregress
from scipy.optimize import curve_fit
from matplotlib import rc_file
from sys import exit
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors

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

fig = plt.figure(figsize=(6.5,5)) 

xmin,xmax = 1.45, 2.45
ymin,ymax = 99,180 

left, right = 0.07, 0.90
bottom, top = 0.08, 0.92
wspace = 0.3

outer = gridspec.GridSpec(1,2,wspace=wspace,left=left,right=right,bottom=bottom,top=top)

fig.text(left +(right-left-wspace/2)/4+0.01,0.99, r"Hydrogen bonding from solvent",ha='center',va='top') 
fig.text(right-(right-left-wspace/2)/4-0.01,0.99, r"Hydrogen bonding from protein",ha='center',va='top') 

fig.text(left +(right-left-wspace/2)/4,0.01, r"$d_{\rm{NH}}$ ($\rm{\AA}$)",ha='center',va='bottom') 
fig.text(right-(right-left-wspace/2)/4,0.01, r"$d_{\rm{NH}}$ ($\rm{\AA}$)",ha='center',va='bottom') 
fig.text(0.01,bottom+(top-bottom)/2         , r"$\theta_1$ (deg)",ha='left',va='center',rotation='vertical') 
fig.text(0.55-wspace/4,bottom+(top-bottom)/2, r"$\theta_1$ (deg)",ha='left',va='center',rotation='vertical') 

fig.text(0.01,0.99,r"\textsf{A}",va='top',ha='left',fontsize=12)
fig.text(0.47,0.99,r"\textsf{B}",va='top',ha='left',fontsize=12)

for subplot in range(2) : 
    #axOuter = fig.add_subplot(outer[subplot]) 
    inner = gridspec.GridSpecFromSubplotSpec(figRows,figCols,subplot_spec=outer[subplot],wspace=0.1,hspace=0.35) 


    for index,molec in enumerate(molecList) : 
        if subplot == 0 : 
            file="SNase_%s/hbond/geometry.xvg"%molec
        elif subplot == 1 : 
            file="SNase_%s/hbond/nw_geometry.xvg"%molec
        else : 
            print "How did you get here? " 
            exit() 

        ax = plt.Subplot(fig,inner[index%figRows,index/figRows]) 
        fig.add_subplot(ax,sharex='all',sharey='all') 
#        ax.set_xlim([xmin,xmax]) 
#        ax.set_ylim([ymin,179]) 

        if not index%figRows + 1 == figRows :
            ax.tick_params(axis="x",labelbottom=False)  
        if index/figRows + 1 == figCols :
            ax.tick_params(axis="y",labelleft=False)  
        ax.text(0.50,1.00,"%s"%molec,va='bottom',ha='center',color=colorDict[molec],transform=ax.transAxes)

#        continue 

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
    
        xbins, ybins = np.arange(xmin, xmax, 0.010), np.arange(ymin,ymax,1.0)
        z, x, y = np.histogram2d(data[:,2],data[:,3],[xbins,ybins])

        im = ax.pcolor(x,y,z.T, cmap='plasma', norm=colors.LogNorm(vmin=1, vmax=10))



#        x = data[:,2] ; y = data[:,3]
#        counts,  xbins,  ybins  = np.histogram2d(x, y,range=((xmin,xmax),(ymin,ymax))) #,bins=(32,32))
#
#        ax.contour(counts.transpose(),extent=[xmin,xmax,ymin,ymax],  
#           linewidths=1,colors='black',  
#           linestyles='solid',levels=np.geomspace(1,10000,10))

    

#cbar_ax = fig.add_axes([0.90, 0.15, 0.015, 0.7]) 
cbar_ax = fig.add_axes([(1 - right)/2 + right - 0.02, bottom+(top-bottom-0.7)/2, 0.015, 0.7]) 
fig.text(0.99,bottom+(top-bottom-0.7)/2 + 0.7 / 2, r"Counts per bin", ha='right', va='center',rotation='vertical') 
cbar = fig.colorbar(im, cax=cbar_ax,ticks=range(1,11) ) 
cbar.set_clim(1,10) 
cbar_ax.set_yticklabels([1,2,3,4,5,6,7,8,9,r'$>$10']) 
#cbar_ax.set_ylim([0,10]) 

fig.savefig('figures/Geometry_combined.png',format='png') 

