import glob 
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.colors import LogNorm
import os
import matplotlib as mpl 
from scipy.stats import linregress
from scipy.optimize import curve_fit
from matplotlib import rc_file
from matplotlib.lines import Line2D
from sys import exit

rcFile = 'rc_files/paper.rc'  
absdata = 'Exp_data/abs_data2.dat' 

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

figCols=3
figRows=4

rc_file(rcFile) 

if not os.path.isdir('figures') : 
    os.mkdir('figures') 


#datafiles = glob.glob('B_State/*/fit_hbond_with_ca/dist.poly') 
#
###Theta 1 analyis

left, right = 0.15, 0.80 
bottom, top = 0.15, 0.95

f2, ax2 = plt.subplots(1,1,figsize=(3.25,2)) #,figsize=(4.3,3) )
f2.subplots_adjust(left=left,bottom=bottom,right=right,top=top) 
f2.text(left+(right-left)/2,0.01, r"$\theta_1$ (deg)", ha='center', va='bottom') 
f2.text(0.01,bottom+(top-bottom)/2, r"$\tilde{\nu}$ (cm$^{-1}$)", ha='left', va='center',rotation='vertical') 

r = 2.45 

peakDict, fwhmDict = {}, {}
with open(absdata) as f : 
    for line in f.readlines() :
        if line.startswith('#') : continue
        key,peak,peakError,fwhm,fwhmError = line.split()
        if fwhm == "nan" : continue
#        value = float(value) * -1
        fwhmDict[key] = [float(fwhm),float(fwhmError) ]
        peakDict[key] = [float(peak),float(peakError) ]

accumAvgTheta = [] 
accumAbsMax = [] 
accumNumHBonds = [] 

legend = [] 

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

    numHBonds = 0 
    if solvent : 
        numHBonds += len(dataWgeom[:,0]) 

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
        numHBonds += len(dataPgeom[:,0]) 

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

    #if not molec == "A90X" : 
    if True : #numHBonds > 1500 : 
        try : 
            ax2.scatter(avgAngleTot,peakDict[molec][0], marker='o',color=colorDict[molec],edgecolor='none',s=numHBonds/50,zorder=10) 
            #ax2.scatter(peakDict[molec],avgAngleTot, marker='o',label=molec,color=colorDict[molec],edgecolor='none',s=numHBonds/50,zorder=10)  
            accumAbsMax.append(peakDict[molec][0]) 
            accumAvgTheta.append(avgAngleTot) 
            accumNumHBonds.append(numHBonds)
            print molec, avgAngleTot

            legend.append(Line2D([0], [0], marker='o', color='none', label=molec,markerfacecolor=colorDict[molec], markersize=5))
        except KeyError : 
            print "No key found for %s"%molec
            continue 
    else : 
            ax2.scatter(avgAngleTot,peakDict[molec][0], marker='o',label=molec,color=colorDict[molec],edgecolor='none',s=numHBonds/50,zorder=10,alpha=0.25)  
            #ax2.scatter(peakDict[molec],avgAngleTot, marker='o',label=molec,color=colorDict[molec],edgecolor='none',s=numHBonds/50,zorder=10,alpha=0.25)  


    #ax2.axhline(2227.5,color='k',linestyle='--') 
    #ax2.axhline(2235.9,color='#66CDFF',linestyle='--') 

ax2.set_xlim(115,160) 

y,x,weights = accumAbsMax, accumAvgTheta, accumNumHBonds
x = np.array(x) 
y = np.array(y) 
weights = np.array(weights) 
#weights = np.ones_like(weights) 
#sigma = 1/weights * 1000

def f(x, m, b): 
    return m*x+b 
def rms(y, yfit):
    return np.sqrt(np.sum((y-yfit)**2))

# Unweighted fit
#p0 = 10, 2160
#p0 = 4.23330709e+00,  -9.00912243e+03
#popt, pcov = curve_fit(f, x, y, p0,sigma=sigma,absolute_sigma=False)
#popt, pcov = curve_fit(f, x, y, p0) 

popt,pcov = np.polyfit(x,y,1,cov=True,w=weights) 

yfit = f(x, *popt)

residuals = y - yfit
ss_res = np.sum((residuals**2)*weights) 
ss_tot = np.sum(((y-np.mean(y))**2)*weights)

r_squared = 1 - (ss_res / ss_tot)
print "R^2 = %0.3f"%((r_squared)) 

print('Unweighted fit parameters:', popt)
print('Covariance matrix:'); print(pcov)
print('rms error in fit:', rms(y, yfit))
print()



#slope,intercept,r_value,p_value,std_error = linregress(accumAvgTheta, accumAbsMax)
slope, intercept = popt
r_value = np.sqrt(r_squared) 
#print "r = %f, p = %f"%(r_value,p_value) 

xs = np.linspace(min(x), max(x) ) 
ys = slope * xs + intercept
ax2.plot(xs, ys, color='k')

ax2.text(0.05,0.90,r"$r$ = %0.2f"%r_value,transform=ax2.transAxes) 
ax2.legend(handles=legend,loc=(0.95,0.10)) 

f2.savefig('figures/weighted_abs_max_vs_max_theta.png',format='png') 
plt.close() 
