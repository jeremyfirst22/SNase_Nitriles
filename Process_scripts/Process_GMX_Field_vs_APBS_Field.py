import numpy as np 
import matplotlib.pyplot as plt 
from scipy.stats import linregress
import sys

nameToColorKeys = {
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
"V23X",
"L25X",
"L38X",
"A58X",
"T62X",
"V66X",
"A90X",
"I92X",
"A109X",
#"V104X",
"N118X"
]

apbs, gmx = [],[]
for molec in molecList : 
    apbsField = 'SNase_%s/APBS_fixed/coloumb_field.out'%molec
    gmxField = 'SNase_%s/force_calc/SNase_%s.protein_field.projected.xvg'%(molec,molec)

    try : 
        apbsData = np.genfromtxt(apbsField)
        gmxData = np.genfromtxt(gmxField) 
    except IOError : 
        continue 
    
    plt.scatter(np.average(apbsData),np.average(gmxData),color=nameToColorKeys[molec]) 
    apbs.append(np.average(apbsData)) 
    gmx.append(np.average(gmxData)) 

slope, intercept, r_value, p_value, std_error = linregress(apbs,gmx) 
xs = np.linspace(np.min(apbs),np.max(apbs),100)
print "%0.2f\t%0.3f"%(r_value,p_value)
plt.plot(xs,slope*xs+intercept,'k--',label="r = %0.3f"%r_value)
plt.xlabel('APBS Field') 
plt.ylabel('GMX Field') 
plt.legend() 

plt.savefig('figures/gmx_pcf_vs_apbs_pcf.png',format='png') 
plt.close()

apbs, gmx = [],[]
for molec in molecList : 
    apbsField = 'SNase_%s/APBS_fixed/rxn_field.out'%molec
    gmxField = 'SNase_%s/force_calc_ca/SNase_%s.solvent_rxn_field.projected.xvg'%(molec,molec)

    try : 
        apbsData = np.genfromtxt(apbsField)
        gmxData = np.genfromtxt(gmxField) 
    except IOError : 
        continue 


    apbsData = np.genfromtxt(apbsField)
    gmxData = np.genfromtxt(gmxField) 
    
    plt.scatter(np.average(apbsData),np.average(gmxData),color=nameToColorKeys[molec]) 
    apbs.append(np.average(apbsData)) 
    gmx.append(np.average(gmxData)) 

slope, intercept, r_value, p_value, std_error = linregress(apbs,gmx) 
xs = np.linspace(np.min(apbs),np.max(apbs),100)
print "%0.2f\t%0.3f"%(r_value,p_value)
plt.plot(xs,slope*xs+intercept,'k--',label="r = %0.3f"%r_value)
plt.xlabel('APBS Field') 
plt.ylabel('GMX Field') 
plt.legend() 

plt.savefig('figures/gmx_srf_vs_apbs_srf.png',format='png') 
plt.close()
