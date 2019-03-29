import glob 
import numpy as np 
import matplotlib.pyplot as plt 
import os
from os import sys
from matplotlib.colors import LogNorm
import matplotlib.lines as mlines 
from matplotlib import rc_file
from scipy.stats import linregress
from scipy.optimize import curve_fit

figCols=1
figRows=1

ymin,ymax = 2156,2166

colorDict = {
        "V23X":'#007F00',
        "L25X":'#4DBEEE',
        "L38X":'#142B8C',
        "A58X":'#BF00BF',
        "T62X":'#77AB30',
        "V66X":'#D95219',
        "A90X":'#A2142F',
        "I92X":'#6666FF',
        "A109X":'#ECB120',
        "V104X":'k',
        "N118X":'#7E2F8E'
}

molecList = [
#"V23X",
#"L25X",
"L38X",
"A58X",
"T62X",
#"V66X",
#"A90X",
#"I92X",
"A109X",
#"V104X",
"N118X"   
]

peakDict = {}  
with open('Exp_data/sasa_abs_data.dat') as f : 
    for line in f : 
        if not line.startswith('#') : 
            key = line.split()[0] 
            value = float(line.split()[2]) 
            peakDict[key] = value 


rcFile='rc_files/paper.rc'
rc_file(rcFile) 

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

if not os.path.isdir('data_fitting') : 
    os.mkdir('data_fitting') 

def findContributions(equilTime) : 
    global angDict, distDict, ang2Dict, srfDict, pcfDict, numHBDict, maxHBonds
    simTime = 50 ## ns 
    angDict, distDict, ang2Dict, srfDict, pcfDict, numHBDict = {},{},{},{},{},{}
    maxHBonds = 0 
    for molec in molecList :
        fileSRF = "SNase_%s/force_calc/SNase_%s.solvent_rxn_field.projected.xvg"%(molec,molec) 
        filePCF = "SNase_%s/force_calc/SNase_%s.protein_field.projected.xvg"%(molec,molec) 
    
        dataSRF = np.genfromtxt(fileSRF) 
        dataPCF = np.genfromtxt(filePCF) 
    
        dataSRF = dataSRF[int(float(equilTime)/simTime *len(dataSRF)):] 
        dataPCF = dataPCF[int(float(equilTime)/simTime *len(dataPCF)):] 

        fileSRF="data_fitting/%s_srf"%molec
        np.savetxt(fileSRF+".xvg",dataSRF)
        command = "~/normal_distribution/tiltAngle -f %s.xvg -o %s.out -p %s.poly -g %s.gaus -t 25 --overwrite"%(fileSRF,fileSRF,fileSRF,fileSRF)
        os.system(command)
        dataSRF = np.genfromtxt('%s.gaus'%fileSRF) 

        filePCF="data_fitting/%s_pcf"%molec
        np.savetxt(filePCF+".xvg",dataPCF)
        command = "~/normal_distribution/tiltAngle -f %s.xvg -o %s.out -p %s.poly -g %s.gaus -t 25 --overwrite"%(filePCF,filePCF,filePCF,filePCF)
        os.system(command)
        dataPCF = np.genfromtxt('%s.gaus'%filePCF) 

        srfDict[molec] = np.average(dataSRF[:,0],weights=dataSRF[:,1]) * (2.57 * -0.67) 
        pcfDict[molec] = np.average(dataSRF[:,0],weights=dataSRF[:,1]) * (2.57 * -0.67) 
    
        fileSolv = "SNase_%s/hbond/geometry.xvg" %(molec) 
        fileProt = "SNase_%s/hbond/nw_geometry.xvg" %(molec) 
    
        dataSolv = np.genfromtxt(fileSolv,skip_header=23) 
        dataProt = np.genfromtxt(fileProt,skip_header=23) 
    
        if len(dataSolv) > 0 : 
            equilFrames = 0 
            for i in range(len(dataSolv)) : 
                if dataSolv[i,0] * 4 /1000 < equilTime : equilFrames += 1 
                else :  break 
            dataSolv = dataSolv[equilFrames:] 
    
        if len(dataProt) > 0 : 
            equilFrames = 0 
            for i in range(len(dataProt)) : 
                if dataProt[i,0] * 4 / 1000 < equilTime : equilFrames += 1 
                else :  break 
            dataProt = dataProt[equilFrames:] 
    
        hbonds = True 
        if len(dataProt) > 1 and len(dataSolv) > 1 : 
            dataTot = np.vstack((dataSolv,dataProt)) 
        elif len(dataProt) > 1 : 
            dataTot = dataProt
        elif len(dataSolv) > 1 : 
            dataTot = dataSolv
        else : 
            hbonds = False
    
        if hbonds : 
            numHBonds = len(dataProt) + len(dataSolv) 
            if numHBonds > maxHBonds : maxHBonds = numHBonds 
            numHBDict[molec] = numHBonds
    
            data = dataTot[:,2]
            distDict[molec] = np.average(data) 
    
            data = dataTot[:,3]
            file="data_fitting/%s_theta1"%molec
            np.savetxt(file+".xvg",data)
            command = "~/normal_distribution/tiltAngle -f %s.xvg -o %s.out -p %s.poly -g %s.gaus -t 25 --overwrite"%(file,file,file,file)
            os.system(command)
            dataPoly = np.genfromtxt('%s.poly'%file) 

            angles = dataPoly[:,0]
            probs = dataPoly[:,1]

            probs *= numHBonds ##un-normalize

            volumes = np.zeros_like(angles) 
            for i in range(len(volumes) ) : 
                if not i == len(volumes) - 1 : 
                    binSize = angles[i+1] - angles[i] 
                else : 
                    binSize = angles[i] - angles[i-1] 
                r = 2.45 
                volumes[i] = 2*np.pi * r**3 / 3 * (-np.cos((angles[i]+binSize/2) * np.pi / 180.)  + np.cos((angles[i]-binSize/2)* np.pi / 180.) )
            probs /= volumes
            
            angDict[molec] = np.average(angles, weights=probs) 
            #angDict[molec] = np.average(data) 

            data = dataTot[:,4]
            ang2Dict[molec] = np.average(data) 
        else : 
            ang2Dict[molec] = 180
            angDict[molec] = 180
            distDict[molec] = 180
            numHBDict[molec] =0

def three_param(X,a,b,c,d,e,g,h) : 
    #s = X[3] / float(maxHBonds) 
    s = np.exp(-g*X[3] / float(maxHBonds)) + h
    return s*a*(X[2] - b) + c *X[1] + d * X[0] + e 


def calculateVtheory() : 
    peakAccum, pcfAccum, srfAccum, angAccum, numHBAccum = [],[],[],[],[]
    for molec in molecList : 
         #if molec == "A90X" : continue 
         peakAccum.append(peakDict[molec]) 
         pcfAccum.append(pcfDict[molec]) 
         srfAccum.append(srfDict[molec]) 
         angAccum.append(angDict[molec]) 
         numHBAccum.append(float(numHBDict[molec]) ) 
    
    bounds = ([0,100,0,0,2000,0,0],[1,180,1,100,3000,500,500]) 
    popt, pcov = curve_fit(three_param,(srfAccum,pcfAccum,angAccum,numHBAccum),peakAccum,bounds=bounds,maxfev=2000) 

    theoryValues = [] 
    for molec in molecList : 
        #if molec == "A90X" : continue 
        theoryValues.append(three_param((srfDict[molec],pcfDict[molec],angDict[molec],numHBDict[molec]),*popt)) 
    return popt,peakAccum,theoryValues


############################################################################
### Main
############################################################################
for equilTime in np.arange(0,50,5) :
    fig, ax = plt.subplots(1,1,figsize=(4.3,3) )
    fig.subplots_adjust(left=0.20,bottom=0.13,right=0.95,top=0.95)
    fig.text(0.5,0.04, "Field ()", ha='center', va='center') 
    fig.text(0.08,0.5, r"Occurances", ha='center', va='center',rotation='vertical') 
    
    findContributions(equilTime)
    popt, peakAccum,theoryValues = calculateVtheory() 
    
    for i in popt : 
        print "%0.2f"%i,
    print "\n",
    
    #molecList.remove("A90X") 
    for i,molec in enumerate(molecList) : 
        ax.scatter(theoryValues[i], peakDict[molec],color=colorDict[molec]) 
    #molecList.append('A90X') 
    
    slope,intercept,r_value,p_value,std_error = linregress(theoryValues, peakAccum) 
    print "r = %f, p = %f" %(r_value,p_value) 
    
    xs = np.linspace(np.min(theoryValues), np.max(theoryValues), 100)
    ys = slope * xs + intercept
    ax.plot(xs, ys, label="r = %0.3f"%r_value,color='k')
    ax.legend(loc=2) 
    
    plt.savefig('figures/curve_fit_%i.png'%equilTime,format='png') 
